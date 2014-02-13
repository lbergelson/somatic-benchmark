package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import java.io.File
import java.util.Calendar
import java.text.SimpleDateFormat
import org.broadinstitute.sting.queue.util.Logging
import scala.io.Source


sealed abstract class EvaluationGroup
case object HCC1143 extends EvaluationGroup
case object HCC1954 extends EvaluationGroup
case object NormalNormal extends EvaluationGroup




case class TumorNormalPair(tumor: File, normal: File, individual: String, cellLine: EvaluationGroup, reference: File)


class TsvReader(file: File){
    val lines = Source.fromFile(file).getLines()
    val header = lines.next().split("\t")
    val columns: Map[String, Int] = header.zipWithIndex.toMap

    def getLines = lines.map(l => new TsvRow(l, columns))

}
class TsvRow(line: String, columns: Map[String, Int]){
    if( columns.size != cells.length) {
        throw new IllegalStateException(s"Number of columns in header doesn't match number of cells: ${columns.size}!=${cells.length}")
    }

    val cells: Array[String] = line.split("\t")
    def lookup(colName: String): String = {
        val columnIndex = columns(colName)
        cells(columnIndex)
    }
}

object TumorNormalPair {

    def readPairsFile(file: File):Set[TumorNormalPair] = {
        val rows = new TsvReader(file).getLines
        val pairs = rows.map{ row =>
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

    def getTimeStamp = {
        val time = Calendar.getInstance.getTime
        val dateFormat = new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss")
        dateFormat.format(time)
    }

    val timeStamp = getTimeStamp

    val identifierString = s"$name-$version-$timeStamp"


}

class RunCallersOnKDB extends QScript with Logging{
    qscript =>

    @Input(fullName="mutation_caller", shortName="c", doc="Script to invoke the mutation caller.")
    var caller: File = _

    @Argument(fullName="name", shortName="n", doc="Name of the mutation caller.")
    var name: String = _

    @Argument(fullName="caller_version", shortName="v", doc="Version of the mutation caller.")
    var version: String = _



    val reference: File = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"

    val LIB_DIR = new File(".")
    val TOOL_DIR = new File(LIB_DIR, "tool-scripts")
    val OUTPUT_DIR = new File("/cga/tcga-gsc/benchmark/data/evaluation/")

    val KDB_ANNOTATE_SCRIPT = new File("kdb_annotate.R")


    val pairs = Seq( new TumorNormalPair("/crsp/qa/picard_aggregation/cancer-exome-val-HCC1143-tn-full-07/12345670177/v1/12345670177.bam",
        "/crsp/qa/picard_aggregation/cancer-exome-val-HCC1143-tn-full-02/12345670180/v1/12345670180.bam",
        "hcc1143-77-80", HCC1143, reference),
        new TumorNormalPair("/crsp/qa/picard_aggregation/cancer-exome-val-HCC1954-tn-full-09/12345670183/v1/12345670183.bam",
            "/crsp/qa/picard_aggregation/cancer-exome-val-HCC1954-tn-full-02/12345670186/v1/12345670186.bam",
            "hcc1954-83-86", HCC1954, reference)
    )


    val cellLines = Map(HCC1143 -> new File("/home/unix/louisb/cga_home/kdb/cga_kdb/mongo/HCC1143_calls.maf"),
                        HCC1954 -> new File("/home/unix/louisb/cga_home/kdb/cga_kdb/mongo/HCC1954_calls.maf"))

    /**
     * Builds the CommandLineFunctions that will be used to run this script and adds them to this.functions directly or using the add() utility method.
     */
    def script(): Unit = {



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

        //gather results
        val aggregators: Map[EvaluationGroup, Collector] = aggregateResults(mutationCallerCmds)
        addAll(aggregators.values)





    }



    def aggregateResults(callers: Seq[MutationCallerInvocation]) = {
        val groups: Map[EvaluationGroup, Seq[MutationCallerInvocation]] =  callers groupBy{
            case caller => caller.pair.cellLine
        }

        val collectors: Map[EvaluationGroup, RunCallersOnKDB.this.type#Collector] = groups.mapValues( callers => new Collector( callers.map( c => c.getEvaluator.getSummaryFile)))
        collectors
    }






    class Collector(@Input val evaluationFiles: Seq[File]) extends InProcessFunction{

        override def run(): Unit = {
            val summaries = evaluationFiles map readEvaluationFile

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
                val snps = resultFromLine("_snps", nameToValue)
                val indels = resultFromLine("_indels", nameToValue)
                new Summary(snps=snps, indels=indels)

            } catch {
                case e: java.io.IOException => logger.error("Couldn't read evaluation file d:", e)
                                               throw e
            }

        }

    }


    class Result(TP: Int,FP: Int, FN: Int, Novel: Int)
    class Summary(snps: Result, indels: Result){

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


        def commandLine = required(tool)+
            required(normal)+
            required(tumor)+
            required(reference)+
            required(outputDir)


        lazy val getEvaluator = {
            pair.cellLine match {
                case HCC1143 => new AnnotateKDB(KDB_ANNOTATE_SCRIPT, calls,  cellLines(HCC1143), outputDir)
                case HCC1954 => new AnnotateKDB(KDB_ANNOTATE_SCRIPT, calls,  cellLines(HCC1954), outputDir)

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


