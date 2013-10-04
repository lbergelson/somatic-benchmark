package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import org.broadinstitute.sting.utils.exceptions.UserException.CouldNotReadInputFile
import org.broadinstitute.sting.queue.library.cga.benchmark.{AbbreviatedFile, AnnotatedBamFile }

class RunBenchmark extends QScript {
  qscript =>

//  @Argument(doc="If this is set than only 1 set of files will be run instead of the complete array.", required=false)
//  var is_test = false

  @Argument(fullName="tool", shortName="t", doc="The name of a tool to run.  A matching script named run<Tool>.sh must be placed in the tool-scripts directory.")
  var tool_names: List[String]  = Nil

  @Argument(fullName="no_false_positives", shortName="nofp", doc="Run false positive analysis.", required=false)
  var no_false_positives: Boolean = false

  @Argument(fullName="no_false_negatives", shortName="nofn", doc="Run false negative analysis.", required=false)
  var no_false_negatives: Boolean = false


  @Argument(fullName="list_tools", shortName="list", doc="List the available tool scripts and exit", required=false)
  var list_tools: Boolean = false

  @Argument(fullName="bamtypes", doc="The bamtypes file.", required=false)
  var bamTypesFile : File = new File("bams.bamtype")

  val LIB_DIR = new File(".")

  val GERMLINE_NAME_TEMPLATE = "%s.bam"
  val GERMLINE_MIX_DIR = new File(LIB_DIR, "data_1g_wgs")


  val referenceFile : File = new File("/home/unix/louisb/cga_home/reference/human_g1k_v37_decoy.fasta")

  lazy val (normals,tumors,spikedTumors) = AnnotatedBamFile.splitByType(AnnotatedBamFile.readBamTypesFile(bamTypesFile))
  lazy val spikedNormal = normals.sortBy(_.length).last


  val SPIKE_DIR = new File(LIB_DIR, "fn_data")

  val TOOL_DIR = new File(LIB_DIR, "tool-scripts")


  def script() {

    if(list_tools) listToolsAndExit()

    val tools = getTools(tool_names)

    if (!no_false_positives) {
        val (_, fpCmds) = getFalsePositiveCommands(tools).unzip
        fpCmds.foreach(add(_))
    }

    if (!no_false_negatives) {
        val (_, spikeCmds ) = getSpikedCommands(tools).unzip
        spikeCmds.foreach(add(_))
    }
  }

  def listToolsAndExit() {
    val files = TOOL_DIR.listFiles()
    val tools = files.filter(file => file.getName.startsWith("run") && file.getName.endsWith(".sh"))
    val names = tools.map(_.getName.drop(3).dropRight(3))

    logger.info("======= Tools ========")
    names.foreach( file => logger.info( "=  %s".format(file) ) )
    logger.info("======================")

    System.exit(0)
  }

  def getTools(names: List[String]):List[AbbreviatedFile] = {

      names.map{name =>
          val toolFile = new File(TOOL_DIR, "run%s.sh".format(name))
          if( ! toolFile.exists() ) throw new CouldNotReadInputFile(toolFile, " file does not exist.")
          new AbbreviatedFile(toolFile, name)
      }
  }




         //invokes <tool> with parameters <normal><tumor><reference><outputDir>
  class ToolInvocation extends  CommandLineFunction with RetryMemoryLimit{
            @Input(doc="The script to run")
            var tool: File = _

            @Input(doc="normal sample bam")
            var normal: File = _

            @Input(doc="tumor sample bam")
            var tumor: File = _

            @Input(doc="reference fasta")
            var reference: File = _

            @Output(doc="output directory")
            var outputDir: File =_

            def commandLine = required(tool)+
                              required(normal)+
                              required(tumor)+
                              required(qscript.referenceFile)+
                              required(outputDir)
  }

  object ToolInvocation {
    def apply(tool:File, normal:File, tumor:File, reference:File, outputDir:File) = {
        val ti = new ToolInvocation
        ti.tool = tool
        ti.normal = normal
        ti.tumor = tumor
        ti.reference = reference
        ti.outputDir = outputDir
        ti.memoryLimit = 4
        ti
    }
  }

  def getFalsePositiveCommands(tools: Traversable[AbbreviatedFile]) = {
    def getPureFalsePositivePairs(normals: Traversable[AbbreviatedFile], tumors: Traversable[AbbreviatedFile]) = {
        for{
            normal <- normals
            tumor <- tumors
            if normal != tumor
        } yield (normal, tumor)
    }

    val pureGermline = getPureFalsePositivePairs(normals, tumors )
    generateCmds(tools, pureGermline, "germline_mix")
 }

 def getSpikedCommands(tools:List[AbbreviatedFile]) = {
    val spiked = spikedTumors.map(tumor => (spikedNormal, tumor))
    generateCmds(tools, spiked, "spiked")
 }


  def generateCmds(toolsToTest: Traversable[AbbreviatedFile], normalTumorPairs: Traversable[(AbbreviatedFile, AbbreviatedFile)], outputDir: File):Traversable[(File, CommandLineFunction)] = {
    def generateCmd(tool: AbbreviatedFile, normal:AbbreviatedFile, tumor: AbbreviatedFile, outputDir: File): (File, CommandLineFunction) ={
        val individualOutputDir = new File(outputDir, "%s_N%S_T%S".format(tool.abbreviation, normal.abbreviation, tumor.abbreviation))
        (individualOutputDir, ToolInvocation(tool=tool, normal=normal, tumor=tumor, reference=referenceFile, outputDir=individualOutputDir) )
    }

    for{
     tool <- toolsToTest
     (normal, tumor) <- normalTumorPairs
    } yield generateCmd(tool, normal, tumor, outputDir)
  }


}

