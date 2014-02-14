package org.broadinstitute.cga.benchmark.queue


import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import java.io.{FileWriter, BufferedWriter, File}
import java.util.Calendar
import java.text.SimpleDateFormat
import org.broadinstitute.sting.queue.util.Logging
import scala.io.Source
import org.broadinstitute.sting.commandline


sealed abstract class EvaluationGroup
case object HCC1143 extends EvaluationGroup
case object HCC1954 extends EvaluationGroup
case object NormalNormal extends EvaluationGroup




case class TumorNormalPair(tumor: File, normal: File, individual: String, cellLine: EvaluationGroup, reference: File)


class TsvReader(file: File) extends Logging{
    val lines = Source.fromFile(file).getLines()


    val header = lines.next().split("\t")
    val columns: Map[String, Int] = header.zipWithIndex.toMap
    logger.debug(s"Header:${header.toString}")
    logger.debug(s"Columns:${columns.toString}")
    def getLines: Iterator[TsvRow] = lines.map{l => new TsvRow(l, columns)}

}
class TsvRow(line: String, columns: Map[String, Int]) extends Logging{
    val cells: Array[String] = line.split("\t")

    if( columns.size != cells.length) {
        throw new IllegalStateException(s"Number of columns in header doesn't match number of cells: ${columns.size}!=${cells.length}")
    }


    logger.debug(s"New TsvRow:\n line:$line \n cells:$cells")
    def lookup(colName: String): String = {
        val columnIndex = columns(colName)
        cells(columnIndex)
    }
}

object TumorNormalPair extends Logging{

    def readPairsFile(file: File):Set[TumorNormalPair] = {
        val rows = new TsvReader(file).getLines
        val pairs = rows.map{ row =>
            logger.debug(s"loading row: $row")
            new TumorNormalPair(
                tumor = new File(row.lookup("tumor")),
                normal = new File(row.lookup("normal")),
                individual = row.lookup("individual"),
                cellLine = row.lookup("evaluationGroup") match {
                    case "HCC1143" => HCC1143
                    case "HCC1954" => HCC1954
                    case "NormalNormal" => NormalNormal
                    case default => throw new Exception(s"$default is not a valid EvaluationGroup")
                },
                reference = new File(row.lookup("reference"))
                )
        }
		pairs.toSet
    }
}


case class MutationCallerInformation(caller: File,name: String, version: String){

    private def getTimeStamp = {
        val time = Calendar.getInstance.getTime
        val dateFormat = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss")
        dateFormat.format(time)
    }

    val timeStamp = getTimeStamp

    val identifierString = s"$name-$version"


}

class RunCallersOnKDB extends QScript with Logging{
    qscript =>

    @Input(fullName="mutation_caller", shortName="c", doc="Script to invoke the mutation caller.")
    var caller: File = _

    @Argument(fullName="name", shortName="n", doc="Name of the mutation caller.")
    var name: String = _

    @Argument(fullName="caller_version", shortName="v", doc="Version of the mutation caller.")
    var version: String = _


	@Argument(fullName="pairs", shortName="p", doc="File listing pairs to run on", required=false)
	var pairs_file: File = new File("testpairs.txt")

    val reference: File = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"

    val LIB_DIR = new File(".")
    val TOOL_DIR = new File(LIB_DIR, "tool-scripts")
    val OUTPUT_DIR = new File("/cga/tcga-gsc/benchmark/data/evaluation/")

    val KDB_ANNOTATE_SCRIPT = new File("kdb_annotate.R")
    val COUNT_FP_SCRIPT = new File("countFP.R")


    val test_pairs = Seq( new TumorNormalPair("/crsp/qa/picard_aggregation/cancer-exome-val-HCC1143-tn-full-07/12345670177/v1/12345670177.bam",
        "/crsp/qa/picard_aggregation/cancer-exome-val-HCC1143-tn-full-02/12345670180/v1/12345670180.bam",
        "hcc1143-77-80", HCC1143, reference)
         )


    val cellLines = Map(HCC1143 -> new File("/home/unix/louisb/cga_home/kdb/cga_kdb/mongo/HCC1143_calls.maf"),
                        HCC1954 -> new File("/home/unix/louisb/cga_home/kdb/cga_kdb/mongo/HCC1954_calls.maf"))

    /**
     * Builds the CommandLineFunctions that will be used to run this script and adds them to this.functions directly or using the add() utility method.
     */
    def script(): Unit = {
		val pairs = TumorNormalPair.readPairsFile(pairs_file)
        logger.info(s"Loaded ${pairs.size} pairs.")
        /*
         * 1 run caller on all cell line data sets
         * 2 evaluate each result against the appropriate kdb entry
         * 3 aggregate like entries
         *
         * 4 run caller on micro bams
         * evaluate results
         *
         * 5 run caller on hapmap mixing experiments, evaluate results
         *
         *
         * 4 store values in table ( or mongo...)
         * 5 update website with new values
         *
         */

        //Run caller on all Tumor-Normal Pairs
        val mutCaller = new MutationCallerInformation(caller, name, version)
        val mutationCallerCmds = pairs map ( pair => new MutationCallerInvocation(mutCaller, pair, OUTPUT_DIR) )
        addAll( mutationCallerCmds )

        //Calculate true positives / false positives for each pair
        val evaluators = mutationCallerCmds map(m => m.getEvaluator)
        addAll(evaluators)


        //print summary
        val evalFiles: Seq[File] = evaluators.map(e => e.summary).toSeq
        logger.info(s"Evaluations produced ${evalFiles.length} files.")
        logger.info(s"Evaluation files are ${evalFiles.map(f => f.toString)}")
        val writeResults = new WriteOutResults(evalFiles, mutationCallerCmds, new File("final.results.txt"))
        add(writeResults)

    }








    class WriteOutResults(evaluationSummaries: Seq[File], mutCallers: Traversable[MutationCallerInvocation], @Output resultsFile: File) extends InProcessFunction {
		@Input
		val inputFiles: List[File] = evaluationSummaries.toList

        def aggregateResults(callers: Traversable[MutationCallerInvocation]) = {
            val groups: Map[EvaluationGroup, Traversable[MutationCallerInvocation]] =  callers groupBy{
                case mutCaller => mutCaller.pair.cellLine
            }

            val collectors: Traversable[Collector] = groups.map{ case (evalGroup, groupedCallers) => new Collector(evalGroup, groupedCallers.toSeq)}
            collectors
		} 

        override def run(): Unit =  {
            val aggregators = aggregateResults(mutCallers)

            val bw = new BufferedWriter(new FileWriter(resultsFile) )


            bw.write(s"${Summary.headerString}\tEvaluationGroup\tCaller\tVersion\tTime\n")

            aggregators.foreach( (collector) => bw.write(collector.resultsToString) )

            bw.close()

        }
        
    }



    class Collector( val evalGroup: EvaluationGroup, val mutCallerGroup:Seq[MutationCallerInvocation]){

        val evaluationFiles: Seq[File] = mutCallerGroup.map(mg => mg.getEvaluator.getSummaryFile).toSeq

        val collectionResults:Seq[Summary] = evaluationFiles map readEvaluationFile

        val resultsToString: String = {
            val resultPairs = (collectionResults, mutCallerGroup).zipped
            val lines = resultPairs.map{ case(summary, mutCaller) =>
                val callerInfo = mutCaller.caller
                //s"${Summary.headerString}\tEvaluationGroup\tCaller\tVersion\tTime\n")
                s"$summary\t$evalGroup\t${callerInfo.name}\t${callerInfo.version}\t${callerInfo.timeStamp}\n"
            }
            lines.mkString("")

        }

        def readEvaluationFile(evalFile: File): Summary = {
            import scala.io.Source
            try{
                val lines = Source.fromFile(evalFile).getLines()
                val header = lines.next().split("\t")

                val countAsStrings: Seq[String] = lines.next().split("\t").toSeq
                val counts = countAsStrings.map(s => s.toInt)

                val nameToValue = (header zip counts).toMap
                def resultFromLine(postfix: String, m:Map[String,Int]):Result ={
                    new Result(FP=m("fp"+postfix),
                               FN=m("fn"+postfix),
                               Novel=m("novel"+postfix),
                               TP= m("tp"+postfix))
                }

                val snps = resultFromLine("_snp", nameToValue)
                val indels = resultFromLine("_indel", nameToValue)
                new Summary(snps=snps, indels=indels)

            } catch {
                case e: java.io.IOException => logger.error("Couldn't read evaluation file:", e)
                                               throw e
            }
        }

    }


    class Result(TP: Int,FP: Int, FN: Int, Novel: Int){
        override def toString = s"$TP\t$FP\t$FN\t$Novel"
    }
    object Result{
        def headerString(postfix: String) = s"tp_$postfix\tfp_$postfix\tfn_$postfix\tnovel_$postfix"
    }
    class Summary(snps: Result, indels: Result){
        override def toString = s"${snps.toString()}\t${indels.toString()}"

    }
    object Summary{
        val headerString = s"${Result.headerString("snp")}\t${Result.headerString("indel")}"

    }



    class MutationCallerInvocation(val caller: MutationCallerInformation, val pair: TumorNormalPair, basedir: File) extends  CommandLineFunction with RetryMemoryLimit{
        @Input(doc="The script to run")
        val tool: File = caller.caller

        @Input(doc="normal sample bam")
        val normal: File = pair.normal

        @Input(doc="tumor sample bam")
        val tumor: File = pair.tumor

        @Input(doc="reference fasta")
        val reference: File = pair.reference

        @Output(doc="output directory")
        val outputDir: File = new File(basedir, caller.identifierString+pair.hashCode())

        @Output(doc="maf file of all calls, set to be outputDir/final_calls.maf")
        val calls: File = new File(outputDir, "final_calls.maf")

        memoryLimit = 4

        def commandLine = required(tool)+
            required(normal)+
            required(tumor)+
            required(reference)+
            required(outputDir)


        lazy val getEvaluator = {
            pair.cellLine match {
                case HCC1143 => new AnnotateKDB(KDB_ANNOTATE_SCRIPT, calls,  cellLines(HCC1143), outputDir)
                case HCC1954 => new AnnotateKDB(KDB_ANNOTATE_SCRIPT, calls,  cellLines(HCC1954), outputDir)
                case NormalNormal => new CountFP(COUNT_FP_SCRIPT, calls, outputDir)
            }
        }
    }




    class RscriptCommandLineFunction(@Input val script:File ) extends CommandLineFunction {

        @Argument(doc="List of commandline arguments")
        var  args: List[String] = Nil

        def commandLine: String = required("Rscript")  +required(script)+ repeat(args)
    }

    trait Evaluator extends CommandLineFunction{
        def getSummaryFile: File
    }


    class CountFP(script:File, @Input nnMaf: File, outputDir: File) extends RscriptCommandLineFunction(script) with Evaluator{
        @Output(doc="Annotation Summary")
        val summary: File = new File(outputDir,"nn.summary_kdb.txt")

        args= List(nnMaf, summary)

        override def getSummaryFile: File = summary

    }

    /**
     * run the kdb_annotate.R script to generate a summary file of fales positives / negatives
     * @param script  location of the R script
     * @param mafToAnnotate the maf to annotate
     * @param kdbMaf the maf to annotate based on
     * @param outputDir  the directory to place the outputs in
     */
    class AnnotateKDB( script: File,
                      @Input val mafToAnnotate: File,
                      @Input val kdbMaf: File,
                      outputDir: File) extends RscriptCommandLineFunction(script) with Evaluator {

        private val outputPrefix = outputDir + "/annotated"

        @Output(doc="Annotation Summary")
        val summary: File = new File(outputDir,"annotated.summary_kdb.txt")

        args = List(mafToAnnotate, kdbMaf, outputPrefix, "WEX")

        override def getSummaryFile: File = summary
    }


}


