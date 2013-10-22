package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import java.io.{FileNotFoundException, IOException, PrintWriter}
import org.broadinstitute.sting.queue.util.Logging
import scala.io.Source
import org.apache.commons.io.IOUtils
import org.broadinstitute.sting.queue.library.cga.benchmark.AnnotatedBamFile
import org.broadinstitute.sting.utils.exceptions.UserException

/**
Traverse the output directories of RunBenchmark and gather results.

It assumes that each output is in it's own directory and the final file is called sample.final.indels.vcf"

*/

class GatherResults extends QScript with Logging{
    qscript =>

    @Input(doc = "False positive test root directories.", required = false)
    var false_positive: Seq[File] = List(new File("germline_mix") )

    @Input(doc = "False negative test root directories.", required = false)
    var false_negative: Seq[File] = List(new File("spiked") )

    @Input(doc = "The bamtypes file", required = false)
    var bamtypes: File = "bams.bamtype"

    @Argument(fullName="no_false_positives", shortName="nofp", doc="Run false positive analysis.", required=false)
    var no_false_positives: Boolean = false

    @Argument(fullName="no_false_negatives", shortName="nofn", doc="Run false negative analysis.", required=false)
    var no_false_negatives: Boolean = false

     def script() {
        def runAnalysis(variantType: String) = {

            val resultsFileName = "final.%s.vcf".format(variantType)

            val fpResults = searchForOutputFiles(false_positive, resultsFileName)
            val fnResults = searchForOutputFiles(false_negative, resultsFileName)
            logger.debug("Fp results:" + fpResults)
            logger.debug("Fn results:" + fnResults)

            if (fpResults.isEmpty) logger.info(s"No false Positive $variantType result file detected")
            if (fnResults.isEmpty) logger.info(s"No false Negative $variantType result file detected")

            if(fpResults.nonEmpty && !no_false_positives) analyzePositives(fpResults, variantType)
            if(fnResults.nonEmpty && !no_false_negatives) analyzeNegatives(fnResults, variantType)



            val makeGraphs = new RscriptCommandLineFunction
            makeGraphs.script = "make_graphs.R"
            makeGraphs.args = List("graphs-%s", "falsePositiveCounts-%s.tsv", "falseNegativeCounts-%s.tsv").map( _.format( variantType ))

            add(makeGraphs)
        }

        runAnalysis("snps")
        runAnalysis("indels")

    }



    def analyzePositives(files: Seq[File], variantType: String) = {
        val counter = new countFalsePositives
        counter.input = files
        counter.output = new File("falsePositiveCounts-%s.tsv".format(variantType))
        add(counter)
    }

    class countFalsePositives extends InProcessFunction{
        @Input(doc="false positive vcfs")
        var input: Seq[File] = Nil

        @Output(doc="false positive result file")
        var output: File = _

        def countOneFile(file: File):Int = {
            import scala.io.Source

            val lines = Source.fromFile(file).getLines()
            lines.foldLeft(0)(countVcfLine)


        }

        def countVcfLine(sum:Int, line: String) ={
            line.startsWith("#") match {
                case true => sum
                case false => sum+1
            }
        }

        def run() {
            val counts = input.par.map(countOneFile).seq
            val metaData = input.map(new DirectoryMetaData(_))
            val results = (metaData, counts).zipped map formatOutputLine

            val header = "Tool\tNormal\tTumor\tFalse_Positives"
            printHeaderAndContents(output, header, results)

        }

        def formatOutputLine(metaData: DirectoryMetaData, count: Int):String = {
            "%s\t%s\t%s\t%s".format(metaData.tool,metaData.normalName, metaData.tumorName, count)
        }
    }

    def analyzeNegatives(files: Seq[File], variantType: String) = {
        val (jobs, diffOuts) = files.map{ file =>
            val metaData = new DirectoryMetaData(file)
            val vcfdiff = new VcfDiff
            vcfdiff.outputPrefix = "%s/vcfout-%s".format(file.getParent, variantType)
            vcfdiff.comparisonVcf = swapExt(metaData.tumor.getParentFile, metaData.tumor, "bam","vcf")
            vcfdiff.vcf = file
            val diffOut = new File( file.getParent, "vcfout-%s.diff.sites_in_files".format(variantType))
            vcfdiff.differenceFile = diffOut
            (vcfdiff, diffOut)
        }.unzip

        jobs.foreach(add(_))

        val stats = new ComputeIndelStats{
            this.sitesFiles = diffOuts.toList
            this.results = new File("falseNegativeCounts-%s.tsv".format(variantType))
        }

        add(stats)

    }


    /**
     * Container for directory level meta data.
     */
    class DirectoryMetaData( inputFile: File) {
        val file = if (inputFile.isDirectory) inputFile else inputFile.getParentFile

        def hasSpikeIn: Boolean = splits.length == 4

        //expecting filename in the format of
        // Indelocator_NDEFGHI_T1234_0.4
        // Indelocator_NDEFGHI_T1234

        private val splits = file.getName.split("_")
        val tool: String = splits(0)

        val normalName: String = splits(1).drop(1)
        val normal: File = DirectoryMetaData.getFromBamtypes(normalName)

        val fraction = if(hasSpikeIn){
            Some(splits(3).toDouble)
        } else
            None


        val tumorName :String = splits(2).drop(1)+(fraction.map("_"+_).getOrElse(""))

        val tumor = DirectoryMetaData.getFromBamtypes(tumorName)

    }

    object DirectoryMetaData{

        def getFromBamtypes(name: String) = {
            bams.getOrElse(name, throw new UserException(s"Tried to lookup ${name} in ${bamtypes}, but it did not exist."))
        }
        //bams loaded from bam maps file
        lazy val bams: Map[String, AnnotatedBamFile] = AnnotatedBamFile.readBamTypesFile(bamtypes).map(bam => bam.abbreviation -> bam).toMap
    }



    class ComputeIndelStats extends InProcessFunction(){
        @Input(doc="all files from vcfdiff")
        var sitesFiles: List[File] = _

        @Output(doc="file to print results in")
        var results: File = _

        def run() {
            val indels = sitesFiles.map(extractResultsFromVcftoolsDiff)

            val header= "Tool\tNormal\tTumor\tFraction\tFP\tFN\tMatched"
            val counts = indels.zip(sitesFiles).map{ pair =>
                val ((first,second), file ) = pair

                val onlyFirst = first.diff(second).size
                val onlySecond = second.diff(first).size
                val matches = first.intersect(second).size
                val metaData = new DirectoryMetaData(file)


                "%s\t%s\t%s\t%s\t%s\t%s\t%s".format(metaData.tool, metaData.normalName, metaData.tumorName,
                                                    metaData.fraction.get, onlyFirst, onlySecond,matches)
            }

            printHeaderAndContents(results, header, counts)

        }
    }

    class VcfDiff extends CommandLineFunction {
        @Input(doc="vcf file to be compared")
        var vcf: File = _

        @Input(doc="vcf file to compare against")
        var comparisonVcf: File =_

        @Argument(doc="prefix string for output files")
        var outputPrefix: String =_

        @Output(doc="output file")
        var differenceFile: File = _

        def commandLine: String = required("vcftools") +
                                  required("--vcf", vcf) +
                                  required("--diff", comparisonVcf) +
                                  required("--out", outputPrefix)
    }

    def searchForOutputFiles(roots: Seq[File], resultsFileName: String) = {
       val finalVcfs = roots.flatMap( searchRootDirectory(_,resultsFileName) )
       finalVcfs
    }

    def searchRootDirectory(dir: File, resultsFileName: String) = {
        val resultDirs: Option[Seq[File]] = Option(dir.listFiles())
        val results = resultDirs.getOrElse(Nil).flatMap(checkForResultFile(_,resultsFileName))
        results
    }

    def checkForResultFile(dir: File, resultsFileName: String = "final.indels.vcf"):Option[File] = {
        val files = dir.listFiles()
        if (files != null) {
            files.find( _.getName == resultsFileName )
        } else {
            None
        }
    }


    def extractResultsFromVcftoolsDiff(vcftoolsOut : File): (List[String], List[String]) = {
        def assignIndels(line: String):(Option[String], Option[String]) = {
                val INDEL_POSITION = 1
                val MATCH_POSITION = 2

                val tokens = line.split('\t')
                val indel = tokens(INDEL_POSITION)

                tokens(MATCH_POSITION) match{
                    case "1" => (Some(indel), None)
                    case "2" => (None, Some(indel))
                    case "B" => (Some(indel), Some(indel))
                    case _ => (None, None)
                }
        }

        val indels: List[(Option[String], Option[String])] = try {
            Source.fromFile(vcftoolsOut).getLines().toList.map(assignIndels)
        } catch {
            case e: IOException =>
            logger.error(e.getMessage)
            List((None, None))
        }

        val (first, second) = indels.unzip
        (first.flatten, second.flatten)

    }

    def printHeaderAndContents(outputFile: File, header: String, contents: Traversable[String])={
        val writer = new PrintWriter(outputFile)
        try{
            writer.println(header)
            contents.foreach(writer.println)
        } finally {
            IOUtils.closeQuietly(writer)
        }
    }

    class RscriptCommandLineFunction extends CommandLineFunction {
        @Input(doc="R script to execute")
        var script: File = _

        @Argument(doc="List of commandline arguments")
        var  args: List[String] = Nil

        def commandLine: String = required("Rscript")  +required(script)+ repeat(args)
    }

}
