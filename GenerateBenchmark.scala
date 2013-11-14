package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.Logging
import org.broadinstitute.sting.queue.function.RetryMemoryLimit

import scala.collection.immutable.Map
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.{MergeSamFiles, SortSam}
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel
import org.broadinstitute.variant.variantcontext.VariantContext
import java.io.{IOException, PrintWriter}
import org.apache.commons.io.{FileUtils, IOUtils}
import java.util
import net.sf.samtools.{SAMReadGroupRecord, SAMFileReader}

import scala.collection.JavaConversions._
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.utils.exceptions.UserException.CouldNotReadInputFile

import org.broadinstitute.sting.queue.library.cga.benchmark.BamType._
import org.broadinstitute.sting.queue.library.cga.benchmark.{AnnotatedBamFile, OutputBamTypeFile}

class GenerateBenchmark extends QScript with Logging {
    qscript =>


    val libDir: File = new File(".")

    @Input(fullName="interval_file", doc="File listing intervals to restrict analysis to.", required = false)
    val intervalFile: File = new File(libDir, "benchmark.interval_list")

    @Input(fullName="goldIndels", doc="Variant file with known good indel locations", required=false)
    val indelFile: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Mills_and_1000G_gold_standard.indels.b37.sites.vcf")

    @Input(fullName = "goldSnps", doc="Variant file with known good snp locations", required=false)
    var snpFile: File = new File("/humgen/gsa-hpprojects/GATK/bundle/current/b37/dbsnp_137.b37.vcf")

    @Input(fullName="reference", shortName="R", doc="Reference fasta.", required= false)
    val referenceFile: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

    @Input(fullName="input_bams", shortName="I", doc = "Base bam files")
    var bams: Seq[File] = Nil

    @Input(fullName="spike_contributor_bam", shortName="spike_bam", doc="Bam file to spike variants in from")
    var spikeContributorBAM: File = _

    @Input(doc = "Directory to locate output files in", shortName = "o", required = false)
    var output_dir: File = new File(libDir)

    @Argument(doc = "Run in test mode which produces a limited set of jobs", required = false)
    var is_test: Boolean = false

    @Argument(doc = "Run without generating spiked data", required = false)
    var no_spike: Boolean = false

    @Argument(fullName = "indel_spike_in", shortName="indels",doc = "Perform indel spike in", required = false)
    var indels: Boolean = true

    @Argument(fullName = "point_mutation_spike_in", shortName = "snps", doc = "Perform point mutation spike in", required = false)
    var snps: Boolean = false

    @Argument(doc = "Number of pieces to split bams into, more pieces will take longer but produce more depths of coverage.", required = false)
    var pieces: Int= 6

    @Argument(doc = "Spike directly into the given bams without any mixing steps.",required = false)
    var spike_in_only: Boolean = false

    @Argument(doc = "a name for a bam file.  These will be used if spike_in_only is specified.  The  number of names specified must be 0 or exactly" +
        "match the number of input bams specified", fullName="bam_name", shortName="name", required=false)
    var bam_names: Seq[String] = Nil

    @Argument(doc = "Stops the computation of depth of coverage information", required = false)
    var no_depths: Boolean = false

    var bamNameToFileMap: Map[String, File] = null





    lazy val primaryIndividual = ReadFromBamHeader.getSingleSampleName(bams)

    lazy val spikeInIndividual = ReadFromBamHeader.getSingleSampleName(Seq(spikeContributorBAM))




    lazy val GENOTYPE_MODEL: GenotypeLikelihoodsCalculationModel.Model = {
        import GenotypeLikelihoodsCalculationModel.Model._
        (indels, snps) match {
            case (true, true) => BOTH
            case (false, true) => SNP
            case (true, false) => INDEL
            case _ => BOTH
        }
    }

    def script() = {

        List(intervalFile, referenceFile, indelFile, snpFile, spikeContributorBAM).foreach(ensureFileExists)

        logger.info("Individuals are %s and %s".format(primaryIndividual, spikeInIndividual))

        val mergedBams = if(!spike_in_only) {

            //the last library in the list is considered the "normal"
            val libraries: Seq[String] = ReadFromBamHeader.getLibraries(bams)
            logger.info("Libraries are: "+  libraries.toString())

            //fracture bams
            val fractureOutDir = new File(output_dir, "data_1g_wgs")
            val (splitBams, fractureCmds) = bams.map(FractureBams.makeFractureJobs(_, referenceFile, libraries, pieces, fractureOutDir)).unzip
            fractureCmds.flatten.foreach(add(_))


            //merge bams
            val (mergedBams, mergers) = MergeBams.makeMergeBamsJobs(fractureOutDir, splitBams.flatten, libraries, pieces)
            addAll(mergers)
            mergedBams

        } else {
            //don't perform fracturing, consider the last bam in the list the "normal"
            //if bams have tags use that as the name
            if(bam_names.length == 0) {
                val reverseBams = bams.reverse
                val tumors  =  reverseBams.tail
                val normal = reverseBams.head
                tumors.map(file => new AnnotatedBamFile(file, TUMOR, file.getName)) :+ new AnnotatedBamFile(normal, NORMAL, normal.getName)

            } else {
                if(bam_names.length != bams.length){
                    throw new UserException.BadInput(s"Number of bams specified (${bams.length}) must exactly match the number of names specified (${bam_names.length}");
                }

                val reverseBamsWithNames: Seq[(GenerateBenchmark.this.type#File, String)] = (bams zip bam_names).reverse
                val tumors  =  reverseBamsWithNames.tail
                val normal = reverseBamsWithNames.head
                tumors.map{ info=>val (file, name) = info; new AnnotatedBamFile(file, TUMOR, name)} :+ new AnnotatedBamFile(normal._1, NORMAL, normal._2)

            }
        }



        val spikedBams: Option[Traversable[AnnotatedBamFile]] = if( !no_spike ){
            //make vcfs
            val vcfDir = new File(output_dir, "vcf_data")
            val vcfMakers = MakeVcfs.makeMakeVcfJobs(vcfDir)
            vcfMakers.foreach(add(_))

            //use SomaticSpike to create false negative test data
            val spikeSitesVCF = new File(vcfDir, "%s_Ref_%s_Het.vcf".format(primaryIndividual, spikeInIndividual) )
            val makeFnCommands = new FalseNegativeSim(spikeSitesVCF, spikeContributorBAM)

            val alleleFractions = if (is_test) Set(.8) else Set(0.04, .1, .2, .4, .8)
            val depths = mergedBams.filter( _.typeOfBam == TUMOR)
            val (spikedBams,falseNegativeCmds) = makeFnCommands.makeFnSimCmds(alleleFractions, depths)
            falseNegativeCmds.foreach(add(_))
            Some(spikedBams)
        } else None



        if(!no_depths){
            //compute depth information for generated bams
            val depthFiles = mergedBams.map(file => swapExt(file.getParent,file,"bam","depth"))

            val sounders = (mergedBams,depthFiles).zipped map( (bamFile,depthFile) => new DepthOfCoverage with GeneratorArguments {
                this.input_file:+= bamFile
                this.out = depthFile
                this.omitDepthOutputAtEachBase = true
            })

            val gatherDepths = new GatherDepths
            gatherDepths.depthFiles = depthFiles
            gatherDepths.coverageFile=new File("collectedCoverage.tsv")

            addAll(sounders)
            add(gatherDepths)
        }

        val recordBamTypes = new OutputBamTypeFile
        recordBamTypes.bams = mergedBams ++ spikedBams.getOrElse(Nil)
        recordBamTypes.bamTypeFile = new File("bams.bamtype")
        add(recordBamTypes)

    }


    def ensureFileExists(file: File) {
        if (!file.exists() ) {
            throw new CouldNotReadInputFile(file)
        }
    }


    trait BaseArguments extends CommandLineFunction with RetryMemoryLimit{

    }

    trait GeneratorArguments extends CommandLineGATK with BaseArguments {
        this.reference_sequence = referenceFile
        this.intervals :+= intervalFile
    }

    class FilterByLibrary extends PrintReads with GeneratorArguments {
        @Argument(doc = "library name")
        var library: String = _

        this.memoryLimit = 2
        this.read_filter ++= List("DuplicateRead", "FailsVendorQualityCheck", "UnmappedRead")
        this.simplifyBAM = true

        override def commandLine = super.commandLine + required("-rf", "LibraryRead") + required("--library", library)

    }

    object FractureBams {

        class SplitBam extends CommandLineFunction with BaseArguments {
            @Input(doc = "bam file to copy header from")
            var headerBam: File = _

            @Input(doc = "bam to be split, the one filtered by library")
            var nameSortedBam: File = _

            @Output(doc = "list of output files")
            var outFiles: List[File] = _

            def commandLine = required("%s/splitBam.pl".format(libDir)) +
                required(headerBam) +
                required(nameSortedBam) +
                repeat(outFiles)
        }

        def makeFractureJobs(bam: File, reference: File, libraries: Traversable[String], pieces: Int, outDir: File): (Traversable[File], Traversable[CommandLineFunction]) = {

            def makeSingleFractureJob(libraryName: String): (List[File], List[CommandLineFunction]) = {

                def getSplitSamNames(library: String, pieces: Int): Traversable[String] = {
                    for (i <- 1 to pieces) yield getSplitFileName(library, i, "sam")
                }

                def getCoordinateSortAndConvertToBam(inputSam: File, outputBam: File): CommandLineFunction = {
                    val sort = new SortSam with BaseArguments
                    sort.memoryLimit = 2
                    sort.input :+= inputSam
                    sort.output = outputBam
                    sort.sortOrder = SortOrder.coordinate
                    sort.createIndex = true
                    sort
                }

                val libraryFiltered = new File(outDir, primaryIndividual + ".original.filtered.%s.bam".format(libraryName))
                val filter = new FilterByLibrary {
                    this.memoryLimit = 2
                    this.library = libraryName
                    this.input_file :+= bam
                    this.out = libraryFiltered
                    this.isIntermediate = true
                }

                val sortedBam = new File(outDir, primaryIndividual + ".original.namesorted.%s.bam".format(libraryName))
                val sort = new SortSam with BaseArguments {
                    this.memoryLimit = 16
                    this.maxRecordsInRam = 4000000
                    this.input :+= libraryFiltered
                    this.output = sortedBam
                    this.sortOrder = net.sf.samtools.SAMFileHeader.SortOrder.queryname
                    this.compressionLevel = 1
                    this.isIntermediate=true
                }

                val split = new SplitBam       //Split one bam into multiple numbered sams
                split.headerBam = bam
                split.nameSortedBam = sortedBam

                val splitSams: List[File] = getSplitSamNames(libraryName, pieces).map(new File(outDir, _)).toList
                split.outFiles = splitSams

                val splitBams: List[File] = splitSams.map((sam: File) => swapExt(outDir, sam, "sam", "bam"))

                val converters = (splitSams, splitBams).zipped map {    //into coordinate order and convert to bam format
                    (samFile, outputBam) =>
                        getCoordinateSortAndConvertToBam(samFile, outputBam)
                }

                (splitBams, List(filter, sort, split) ++ converters)
            }

            val (splitBams, cmds) = (for (library <- libraries) yield makeSingleFractureJob(library)).unzip

            (splitBams.flatten, cmds.flatten)
        }
    }

    object MergeBams {
        private class BamFileNameCalculator(val is_test:Boolean, val splitBams: Traversable[File],
                                            val libraries: Seq[String], pieces: Int) {

            val bamNameToFileMap = splitBams.map((bam: File) => (bam.getName, bam)).toMap

            val tumorLibraries = libraries.dropRight(1)
            val normalLibraries = libraries.takeRight(1)

            val normalFiles = calculateFileNames(normalLibraries, is_test)
            val tumorFiles = calculateFileNames(tumorLibraries, is_test)


            /**
             * Returns a list of bams that correspond to the encoded digitString.
             * Each digit in the string maps to specific bam.  Each library is split into multiple files and there is a unique digit
             * assigned to each file.
             * @param digitString
             * @return
             */
            def getBams(digitString: String): List[File] = {
                val bamDigitToNameMap = generateBamMap
                try {
                    digitString.map(digit => bamNameToFileMap(bamDigitToNameMap(digit))).toList
                } catch {
                    case e: Exception =>
                        println(bamNameToFileMap)
                        println(bamDigitToNameMap)
                        throw e
                }
            }



            private def calculateFileNames(libraries: Seq[String], isTest: Boolean) = {
                def allLengthsOfSlice(seq: Seq[Char])={
                    for {
                        length <- 1 to seq.length
                    } yield seq.slice(0, length).mkString("")
                }

                def allLibrarySplitFiles(libraries: Seq[String]):Seq[Char] = {
                    for {
                        library  <- libraries
                        piece <- 1 to pieces
                    } yield calculateDigit(library, piece)
                }


                val files = allLibrarySplitFiles(libraries)
                if(isTest) List(files.mkString("")) else allLengthsOfSlice(files)

            }

            /** Convert the combination of library / string into a unique character
              */
            private def calculateDigit(library: String, piece: Int): Char = {
                val digits = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                val index = libraries.indexOf(library) * pieces + piece - 1
                digits.charAt(index)
            }

            private def generateBamMap: Map[Char, String] = {

                if (pieces * libraries.length > 35) throw new UserException.BadArgumentValue("pieces", "Pieces * number of libraries must be <= 35 due to an implementation detail.  " +
                    "(Pieces:%d, Libraries:%d, P*L=%d.)%n If this is an issue please contact the maintainer.".format(pieces,libraries.length,pieces*libraries.length))

                def mapLibraryPieces[T](f: (String, Int) => T): Seq[T] = {
                    for {
                        library <- libraries
                        piece <- 1 to pieces
                    } yield f(library, piece)
                }

                val mappings = mapLibraryPieces((library, piece) => (calculateDigit(library, piece), getSplitFileName(library, piece, "bam")))

                mappings.toMap

            }
        }

        private val nameTemplate = ".%s.bam"

        def makeMergeBamsJobs(dir: File, splitBams: Traversable[File], libraries: Seq[String], pieces: Int): (Seq[AnnotatedBamFile], Seq[MergeSamFiles]) = {
            val nameCalculator = new BamFileNameCalculator(is_test, splitBams, libraries, pieces )

            def makeJobs(bamNames: Seq[String], bamType: BamType): Seq[(AnnotatedBamFile, MergeSamFiles)] ={
                logger.debug("names="+ bamNames.mkString(",%n".format()))

                bamNames.map {
                    name =>
                        val mergedFile = new AnnotatedBamFile( new File(dir, nameTemplate.format(name)), bamType, name)
                        val inputBams = nameCalculator.getBams(name)
                        val merge = new MergeSamFiles

                        merge.memoryLimit = 2
                        merge.input ++= inputBams
                        merge.output = mergedFile
                        merge.createIndex = true
                        merge.USE_THREADING = true
                        (mergedFile, merge)
                }
            }

            val tumorJobs = makeJobs(nameCalculator.tumorFiles, TUMOR)
            val normalJobs = makeJobs(nameCalculator.normalFiles, NORMAL)
            (tumorJobs++normalJobs).unzip



        }

    }



    class FalseNegativeSim(spikeSitesVCF: File, spikeInBam: File) {
        val spikedOutputDir = new File(qscript.output_dir, "fn_data")

        def makeFnSimCmds(alleleFractions: Traversable[Double], depths: Traversable[AnnotatedBamFile])= {
            val results = for {
                fraction <- alleleFractions
                depth <- depths
            } yield makeMixedBam(fraction, depth)
            val (spikedBams, spikers, alleleCounters) = results.unzip3
            (spikedBams, spikers++alleleCounters)

        }

        private def makeMixedBam(spikeFraction: Double, depth: AnnotatedBamFile): (AnnotatedBamFile, CommandLineFunction, CommandLineFunction) = {
            val tumorBam = depth
            val outBam = new AnnotatedBamFile(new File(spikedOutputDir, deriveBamName(spikeFraction, depth.abbreviation)), SPIKED, deriveBamName(spikeFraction, depth.getName))

            val outIntervals = swapExt(spikedOutputDir, outBam, "bam", "interval_list")
            val outVcf = swapExt(spikedOutputDir, outBam, "bam", "vcf")

            val spiker = new SomaticSpike with GeneratorArguments
            spiker.javaMemoryLimit = 4
            spiker.simulation_fraction = spikeFraction
            spiker.out = outBam
            spiker.input_file :+= tumorBam
            spiker.spiked_intervals_out = outIntervals
            spiker.intervals = List(spikeSitesVCF)
            spiker.input_file :+= new TaggedFile(spikeContributorBAM, "spike")
            spiker.variant = spikeSitesVCF
            spiker.spiked_variants = outVcf

            val alleleCounter = runUgOnSpikedBam(outBam, outVcf, outIntervals)

            (outBam, spiker, alleleCounter)
        }

        private def runUgOnSpikedBam(spikedBam: File, spikedInVariants: File, intervals: File) = {
            val ug = new UnifiedGenotyper with GeneratorArguments
            ug.input_file :+= spikedBam
            ug.alleles = spikedInVariants
            ug.genotype_likelihoods_model = qscript.GENOTYPE_MODEL
            ug.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
            ug.intervals = Seq(intervals)
            ug.out = swapExt(spikedOutputDir, spikedBam, "bam","ug.vcf")
            ug
        }

        private def deriveBamName(alleleFraction: Double, depth: String): String = {
            val bamNameTemplate = "%s_%s_spikein.bam"
            bamNameTemplate.format(depth, alleleFraction)
        }
    }



    def getSplitFileName(library: String, piece: Int, extension: String) = {
        val fileNameTemplate = primaryIndividual + ".split.%s.%03d.%s"
        fileNameTemplate.format(library, piece, extension)
    }

    class GatherDepths extends InProcessFunction {

        @Input(doc="depth files")
        var depthFiles: Seq[File] = Nil

        @Output(doc="coverage summary file")
        var coverageFile: File = _

        def printHeaderAndContents(outputFile: File, header: String, contents: Traversable[String])={
            val writer = new PrintWriter(outputFile)
            try{
                writer.println(header)
                contents.foreach(writer.println)
            } finally {
                IOUtils.closeQuietly(writer)
            }
        }

        def extractAverageCoverage(file:File):String = {
            import scala.collection.JavaConversions._
            val COVERAGE_POSITION = 2
            val lines: util.List[String] = FileUtils.readLines(file)
            val totalLine = lines.find(_.startsWith("Total"))
            val averageCoverage = totalLine.get.split("\t")(COVERAGE_POSITION)
            averageCoverage

        }

        def extractFilename(file:File):String = {
            val splitName = swapExt(file.getParentFile, file, ".depth.sample_summary","").getName.split('.')
            logger.debug("extractFilename: Extracting name from %s: splits as %s".format(file.getName, splitName))
            splitName(splitName.length-1)
        }

        def run() {
            val summaryFiles = depthFiles.map{ file => new File(file.getAbsolutePath+".sample_summary")}
            val filenames = summaryFiles.map(extractFilename)
            val coverage = summaryFiles.map(extractAverageCoverage )

            val header = "File\tCoverage"
            val results = (filenames, coverage).zipped map( (name, coverage) => "%s\t%s".format(name,coverage) )

            printHeaderAndContents(coverageFile, header, results)

        }
    }

    object MakeVcfs {

        def makeMakeVcfJobs(outputDir: File):List[CommandLineFunction] = {
            val genotyperOutputVCF = new File(outputDir, "genotypes_at_known_sites.vcf")

            def makeUnifiedGenotyperJob = {
                val genotyper = new UnifiedGenotyper with GeneratorArguments {
                    this.scatterCount=8
                    this.input_file :+= qscript.spikeContributorBAM
                    this.input_file = input_file ++ qscript.bams
                    this.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
                    this.alleles = new TaggedFile(indelFile, "VCF")
                    this.o = genotyperOutputVCF
                    this.genotype_likelihoods_model = GENOTYPE_MODEL
                }
                genotyper
            }


            def makeSelectVariants(outputVCF: File, selection: Seq[String]) = {
                val selector = new SelectVariants with GeneratorArguments{
                    this.selectType :+= VariantContext.Type.INDEL
                    this.variant = new TaggedFile(genotyperOutputVCF, "VCF")
                    this.select = selection
                    this.out= outputVCF
                }
                selector
            }

            class WriteIntervals extends CommandLineFunction with BaseArguments {
                @Input(doc="vcf file to extract intervals from")
                var inputVCF: File = _

                @Output(doc="interval file")
                var intervals: File = _

                def commandLine: String = required("cat", inputVCF) +
                    required("|", escape = false) +
                    required("grep","-v","#") +
                    required("|", escape = false) +
                    required("awk","{ print $1 \":\" $2 }") +
                    required(">", escape = false) +
                    required(intervals)
            }

            val refAndHetVcf = new File(outputDir, primaryIndividual+"_Ref_"+spikeInIndividual+"_Het.vcf")
            val hetOrHomVcf = new File(outputDir, primaryIndividual+"_Het_Or_HomeNonRef.vcf")

            val genotyper = makeUnifiedGenotyperJob
            val selectFirstRefSecondHet = makeSelectVariants(refAndHetVcf, Seq( "vc.getGenotype(\""+primaryIndividual+"\").isHomRef()"
                                                                               +" && vc.getGenotype(\""+spikeInIndividual+"\").isHet()"
                                                                               +" && vc.getGenotype(\""+spikeInIndividual+"\").getPhredScaledQual() > 50"
                                                                               +" && QUAL > 50") )
            val selectHetOrHomNonRef = makeSelectVariants(hetOrHomVcf, Seq("!vc.getGenotype(\""+primaryIndividual+"\").isHomRef()"
                                                                          +" && vc.getGenotype(\""+primaryIndividual+"\").getPhredScaledQual() > 50"
                                                                          +" && QUAL > 50" ) )



            val writeIntervals = new WriteIntervals{
                this.inputVCF = hetOrHomVcf
                this.intervals = swapExt(outputDir, hetOrHomVcf, "vcf", "intervals")
            }

            List(genotyper, selectFirstRefSecondHet, selectHetOrHomNonRef, writeIntervals)
        }

    }


    /**
     *  Read and validate data from a bam file header.
     */
    object ReadFromBamHeader{

        def getLibraries(bams: Seq[File]): Seq[String] = {
            val libraries = bams.flatMap(getLibraryNames).distinct
            libraries.size match{
                case 0 | 1 => throw new UserException.BadInput("At least 2 libraries are required." )
                case _ => libraries
            }
        }


        private def getLibraryNames(bam: File)= {
            valuesFromReadGroups(bam, record => record.getLibrary, "Couldn't determine Library names.")
        }

        def getSingleSampleName(bams: Seq[File]): String = {
            bams.map(getSampleName).distinct match {
                case names if names.length == 1 => names(0)
                case names => throw new UserException.BadInput("Expected only 1 individual in input bams, but found %d. (%s)"
                    .format(names.length, names.mkString(",")))
            }
        }

        private def getSampleName(bam: File): String = {
                val sample = valuesFromReadGroups(bam, record=> record.getSample, "Couldn't determine Individual name.")
                sample.length match {
                    case 0 => throw new UserException.BadInput("Could not determine any sample name in %s".format(bam.getAbsolutePath))
                    case 1 => sample.head
                    case num => throw new UserException.BadInput("Expecting only a single sample from %s, found %d".format(bam.getAbsolutePath, num))
                }

        }

        private def valuesFromReadGroups[A] (bam: File, getter: (SAMReadGroupRecord => A), messageOnFailure: String)={
            var reader: SAMFileReader = null

            try {
                reader = new SAMFileReader(bam)
                val values = reader.getFileHeader.getReadGroups.map(getter).distinct
                reader.close()
                values

            } catch {
                case e: IOException => throw new CouldNotReadInputFile(bam,messageOnFailure, e)
            } finally {
                IOUtils.closeQuietly(reader)
            }

        }
    }


}
