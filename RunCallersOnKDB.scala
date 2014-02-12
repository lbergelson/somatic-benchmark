package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import java.io.File
import java.util.Calendar
import java.text.SimpleDateFormat


sealed abstract class EvaluationGroup
case object HCC1143 extends EvaluationGroup
case object HCC1954 extends EvaluationGroup
case object NormalNormal extends EvaluationGroup




case class TumorNormalPair(tumor: File, normal: File, individual: String, cellLine: EvaluationGroup, reference: File)


class MutationCallerInformation(val caller: File, val name: String, val version: String){

    def getTimeStamp = {
        val time = Calendar.getInstance.getTime
        val dateFormat = new SimpleDateFormat("yyyy-MM-dd-HH:mm:ss")
        dateFormat.format(time)
    }

    val timeStamp = getTimeStamp

    val identifierString = s"$name-$version-$timeStamp"


}

class RunCallersOnKDB extends QScript{
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

    val KDB_ANNOTATE_SCRIPT = new File(LIB_DIR, "kdb_annotate.R")


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
        val evaluators = mutationCallerCmds map getEvaluator
        addAll(evaluators)





    }

    def getEvaluator(caller: MutationCallerInvocation) = {
        caller.pair.cellLine match {
            case HCC1143 => new AnnotateKDB(KDB_ANNOTATE_SCRIPT, caller.calls,  cellLines(HCC1143), caller.outputDir)
            case HCC1954 => new AnnotateKDB(KDB_ANNOTATE_SCRIPT, caller.calls,  cellLines(HCC1954), caller.outputDir)

        }
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
        val outputDir: File = new File(basedir, caller.identifierString)

        @Output(doc="maf file of all calls, set to be outputDir/final_calls.maf")
        val calls: File = new File(outputDir, "final_calls.maf")


        def commandLine = required(tool)+
            required(normal)+
            required(tumor)+
            required(reference)+
            required(outputDir)


    }




    class RscriptCommandLineFunction(@Input val script:File ) extends CommandLineFunction {

        @Argument(doc="List of commandline arguments")
        var  args: List[String] = Nil

        def commandLine: String = required("Rscript")  +required(script)+ repeat(args)
    }

    class AnnotateKDB( script: File,
                      @Input val mafToAnnotate: File,
                      @Input val kdbMaf: File,
                      outputDir: File) extends RscriptCommandLineFunction(script) {

        val outputPrefix = outputDir + "/annotated"

        @Output(doc="Annotation Summary")
        val summary: File = new File(outputPrefix,".summary_kdb.txt")

        args = List(mafToAnnotate, kdbMaf, outputPrefix, "WEX")

    }
}


