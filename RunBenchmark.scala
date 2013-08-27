package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.io.{IOUtils, FileExtension}
import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import org.broadinstitute.sting.utils.exceptions.UserException.CouldNotReadInputFile
import java.io.{IOException, FileReader}
import scala.io.Source
import org.apache.commons.io.IOUtils

class RunBenchmark extends QScript {
  qscript =>

//  @Argument(doc="If this is set than only 1 set of files will be run instead of the complete array.", required=false)
//  var is_test = false

  @Argument(fullName="tool", shortName="t", doc="The name of a tool to run.  A matching script named run<Tool>.sh must be placed in the tool-scripts directory.")
  var tool_names: List[String]  = Nil

  val GERMLINE_MIX_DIR = new File("data_1g_wgs")
  val SPIKE_DIR = new File("fn_data")


  val bamTypesFile : File = new File("bams.bamtype")
  val referenceFile : File = new File("/home/unix/louisb/cga_home/reference/human_g1k_v37_decoy.fasta")

  lazy val (normals,tumors,spikedTumors) = readBamTypesFile(bamTypesFile)
  lazy val spikedNormal = normals.sortBy(_.length).last




  def script() {


    val tools = getTools(tool_names)
    val (outFPDirs, fpCmds) = getFalsePositiveCommands(tools).unzip
    val (outSpikedDirs, spikeCmds ) = getSpikedCommands(tools).unzip

    (fpCmds ++ spikeCmds).foreach(add(_))
  }

  def getTools(names: List[String]):List[AbrvFile] = {
      val toolDir = "tool-scripts"
      names.map{name =>
          val toolFile = new File(toolDir, "run%s.sh".format(name))
          if( ! toolFile.exists() ) throw new CouldNotReadInputFile(toolFile, " file does not exist.")
          new AbrvFile(toolFile, name)
      }
  }

  class AbrvFile(file: File , val abrv: String) extends File(file) with FileExtension {
    def withPath(path: String) = new AbrvFile(path, abrv)
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

  def getFalsePositiveCommands(tools: Traversable[AbrvFile]) = {
    def getPureFalsePositivePairs(normals: Traversable[AbrvFile], tumors: Traversable[AbrvFile]) = {
        for{
            normal <- normals
            tumor <- tumors
        } yield (normal, tumor)
    }

    val pureGermline = getPureFalsePositivePairs(normals, tumors )
    generateCmds(tools, pureGermline, "germline_mix")
 }

 def getSpikedCommands(tools:List[AbrvFile]) = {
    val spiked = spikedTumors.map(tumor => (spikedNormal, tumor))
    generateCmds(tools, spiked, "spiked")
 }


  def generateCmds(toolsToTest: Traversable[AbrvFile], normalTumorPairs: Traversable[(AbrvFile, AbrvFile)], outputDir: File):Traversable[(File, CommandLineFunction)] = {
    def generateCmd(tool: AbrvFile, normal:AbrvFile, tumor: AbrvFile, outputDir: File): (File, CommandLineFunction) ={
        val individualOutputDir = new File(outputDir, "%s_N%S_T%S".format(tool.abrv, normal.abrv, tumor.abrv))
        (individualOutputDir, ToolInvocation(tool=tool, normal=normal, tumor=tumor, reference=referenceFile, outputDir=individualOutputDir) )
    }

    for{
     tool <- toolsToTest
     (normal, tumor) <- normalTumorPairs
    } yield generateCmd(tool, normal, tumor, outputDir)
  }

  def readBamTypesFile(typesFile: File) = {
    try{

        val source = Source.fromFile(typesFile)

        val lines = source.getLines()
        val splits = lines.map(_.split('\t')).map(splitLine => (splitLine(0), splitLine(1))).toList

        def newGermlineMixFile(name : String)= {
            new AbrvFile(new File(GERMLINE_MIX_DIR, name), name.split('.')(1) )
        }
        def newSpikedFile(name : String) = {
            new AbrvFile(new File(SPIKE_DIR, name), name.split('_')(1))
        }

        val tumor = splits.collect  { case ("TUMOR", name) => newGermlineMixFile(name)  }
        val normal = splits.collect { case ("NORMAL", name) => newGermlineMixFile(name) }
        val spiked = splits.collect { case ("SPIKED", name) => newSpikedFile(name) }

        logger.debug("Normal: "+normal)
        logger.debug("Tumor: "+tumor)
        logger.debug("Spiked: " + spiked)

        (normal, tumor, spiked)
    } catch {
        case e: IOException => throw new CouldNotReadInputFile(typesFile, "Could not read bam types file", e)
    }
  }

}

