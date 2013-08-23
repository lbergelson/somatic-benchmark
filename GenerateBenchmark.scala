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
import net.sf.samtools.SAMFileReader

import scala.collection.JavaConversions._
import org.broadinstitute.sting.utils.exceptions.UserException

class GenerateBenchmark extends QScript with Logging {
    qscript =>

    //TODO implement these as cmdline parameters instead of hard coding them
    val indelFile: File = new File("/humgen/1kg/DCC/ftp/technical/working/20130610_ceu_hc_trio/broad/CEU.wgs.HaplotypeCaller_bi.20130520.snps_indels.high_coverage_pcr_free.genotypes.vcf.gz")
    val referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

    val libDir: File = new File(".")

    @Input(fullName ="input_bams", shortName="I", doc = "Base bam files")
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



    var bamNameToFileMap: Map[String, File] = null

    val PIECES = 6

    //the last library in the list is considered the "normal"
    lazy val libraries = getLibraries(bams)

    lazy val primaryIndividual = bams.map(getIndividualName).distinct match {
        case names if names.length == 1 => names(0)
        case names => throw new UserException.BadInput("Expected only 1 individual in input bams, but found %d. (%s)"
            .format(names.length, names.mkString(",")))
    }
    lazy val spikeInIndividual = getIndividualName(spikeContributorBAM)

    val intervalFile = new File(libDir, "benchmark.interval_list")

    lazy val SAMPLE_NAME_PREFIX = primaryIndividual+".WGS"
    val prefix = "chr1"

    lazy val tumorFiles = allLibrarySplitFiles(tumorLibraries)
    lazy val normalFiles = allLibrarySplitFiles(normalLibraries)

    def script() = {

        logger.info("Libraries are: "+  libraries.toString())
        logger.info("Individuals are %s and %s".format(primaryIndividual, spikeInIndividual))


        //fracture bams
        val fractureOutDir = new File(output_dir, "data_1g_wgs")
        val (splitBams, fractureCmds) = bams.map(FractureBams.makeFractureJobs(_, referenceFile, libraries, PIECES, fractureOutDir)).unzip
        fractureCmds.flatten.foreach(add(_))

        qscript.bamNameToFileMap = splitBams.flatten.map((bam: File) => (bam.getName, bam)).toMap

        if( !no_spike ){
            //make vcfs
            val vcfDir = new File(output_dir, "vcf_data")
            val vcfMakers = MakeVcfs.makeMakeVcfJobs(vcfDir)
            vcfMakers.foreach(add(_))

            //use SomaticSpike to create false negative test data
            val spikeSitesVCF = new File(vcfDir, "%s_Ref_%s_Het.vcf".format(primaryIndividual, spikeInIndividual) )
            val makeFnCommands = new FalseNegativeSim(spikeSitesVCF, spikeContributorBAM)

            val alleleFractions = Set(0.04, .1, .2, .4, .8)
            val depths = if(is_test) List(tumorFiles.mkString("")) else allLengthsOfSlice(tumorFiles)
            val (_,falseNegativeCmds) = makeFnCommands.makeFnSimCmds(alleleFractions, depths)
            falseNegativeCmds.foreach(add(_))
        }

        //merge bams
        val (mergedBams, mergers) = MergeBams.makeMergeBamsJobs(fractureOutDir)
        addAll(mergers)

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

        def makeFractureJobs(bam: File, reference: File, libraries: Traversable[String], pieces: Int, outDir: File) = {

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
        private val outFileNameTemplate = primaryIndividual+".%s.bam"
        private lazy val BAMGROUPS: Seq[String] = if (is_test) {
            List(tumorFiles.mkString(""), normalFiles.mkString(""))
        } else {
            val tumors: Seq[String] = allLengthsOfSlice(tumorFiles)
            val normals: Seq[String] = allLengthsOfSlice(normalFiles)
            tumors++normals
        }


        def makeMergeBamsJobs(dir: File) = {

            logger.debug("BAMGROUPS="+ BAMGROUPS.mkString(",%n".format()))
            BAMGROUPS.map {
                name =>
                    val mergedFile = new File(dir, outFileNameTemplate.format(name))
                    val inputBams = getBams(name)
                    val merge = new MergeSamFiles

                    merge.memoryLimit = 2
                    merge.input ++= inputBams
                    merge.output = mergedFile
                    merge.createIndex = true
                    merge.USE_THREADING = true
                    (mergedFile, merge)
            }.unzip
        }


    }



    class FalseNegativeSim(spikeSitesVCF: File, spikeInBam: File) {
        val spikedOutputDir = new File(output_dir, "fn_data")

        def makeFnSimCmds(alleleFractions: Traversable[Double], depths: Traversable[String])= {
            val pairs = for {
                fraction <- alleleFractions
                depth <- depths
            } yield makeMixedBam(fraction, depth)
            pairs.unzip
        }

        private def makeMixedBam(alleleFraction: Double, depth: String): (File, CommandLineFunction) = {
            val tumorBams = getBams(depth)
            val outBam = new File(spikedOutputDir, deriveBamName(alleleFraction, depth))
            val outIntervals = swapExt(spikedOutputDir, outBam, "bam", "interval_list")
            val outVcf = swapExt(spikedOutputDir, outBam, "bam", "vcf")

            val spike = new SomaticSpike with GeneratorArguments
            spike.javaMemoryLimit = 4
            spike.simulation_fraction = alleleFraction
            spike.out = outBam
            spike.input_file ++= tumorBams
            spike.spiked_intervals_out = outIntervals
            spike.intervals == List(spikeSitesVCF)
            spike.input_file :+= new TaggedFile(spikeContributorBAM, "spike")
            spike.variant = spikeSitesVCF
            spike.spiked_variants = outVcf
            (outBam, spike)
        }

        private def deriveBamName(alleleFraction: Double, depth: String): String = {
            val bamNameTemplate = primaryIndividual+"_%s_"+spikeInIndividual+"_%s_spikein.bam"
            bamNameTemplate.format(depth, alleleFraction)
        }
    }

    def tumorLibraries = {
        libraries.dropRight(1)
    }

    def normalLibraries = {
        libraries.takeRight(1)
    }

    def allLengthsOfSlice(seq: Seq[Char])={
        for {
            length <- 1 to seq.length
        } yield seq.slice(0, length).mkString("")
    }

    def allLibrarySplitFiles(libraries: Seq[String]):Seq[Char] = {
        for {
            library  <- libraries
            piece <- 1 to PIECES
        } yield calculateDigit(library, piece)
    }

    /**
     * Returns a list of bams that correspond to the encoded digitString.
     * Each digit in the string maps to specific bam.  Each library is split into multiple files and there is a unique digit
     * assigned to each file.
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

    /** Convert the combination of library / string into a unique character
      */
    def calculateDigit(library: String, piece: Int): Char = {
        val digits = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        val index = libraries.indexOf(library) * PIECES + piece - 1
        digits.charAt(index)
    }



    def generateBamMap: Map[Char, String] = {

        if (PIECES > 6 || libraries.size > 3) throw new UnsupportedOperationException("Currently only supported for PIECES <= 6 and LIBRARIES.size <= 3")

        def mapLibraryPieces[T](f: (String, Int) => T): Seq[T] = {
            for {
                library <- libraries
                piece <- 1 to PIECES
            } yield f(library, piece)
        }

        val mappings = mapLibraryPieces((library, piece) => (calculateDigit(library, piece), getSplitFileName(library, piece, "bam")))

        mappings.toMap

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
                    this.scatterCount=4
                    this.input_file :+= qscript.spikeContributorBAM
                    this.input_file = qscript.bams
                    this.genotyping_mode = GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
                    this.alleles = new TaggedFile(indelFile, "VCF")
                    this.o = genotyperOutputVCF
                    this.genotype_likelihoods_model = (qscript.indels,qscript.snps) match {
                        case (true, true) => GenotypeLikelihoodsCalculationModel.Model.BOTH
                        case (false, true) => GenotypeLikelihoodsCalculationModel.Model.SNP
                        case (true, false) => GenotypeLikelihoodsCalculationModel.Model.INDEL
                        case _ => GenotypeLikelihoodsCalculationModel.Model.BOTH
                    }
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

    def getLibraries(bams: Seq[File]) = {
        val libraries = bams.flatMap(getLibraryNames).distinct
        libraries.size match{
            case 0 | 1 => throw new UserException.BadInput("At least 2 libraries are required." )
            case _ => libraries
        }
    }


    def getLibraryNames(bam: File)= {
        try{
            val reader: SAMFileReader = new SAMFileReader(bam)
            val libraries = reader.getFileHeader.getReadGroups.map(_.getLibrary).distinct
            libraries
        } catch {
            case e: IOException => throw new UserException.CouldNotReadInputFile(bam, "Couldn't determine library names.", e)
        }
    }


    def getIndividualName(bam: File)= {
        try {
            val reader: SAMFileReader = new SAMFileReader(bam)
            val sample = reader.getFileHeader.getReadGroups.map(_.getSample).distinct

            sample.length match {
                case 0 => throw new UserException.BadInput("Could not determine any sample name in %s".format(bam.getAbsolutePath))
                case 1 => sample.head
                case num => throw new UserException.BadInput("Expecting only a single sample from %s, found %d".format(bam.getAbsolutePath, num))
            }

        } catch {
            case e: IOException => throw new UserException.CouldNotReadInputFile(bam, "Couldn't determine sample name.", e)
        }
    }


}


