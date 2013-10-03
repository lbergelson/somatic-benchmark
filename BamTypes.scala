package org.broadinstitute.cga.benchmark.queue

import org.broadinstitute.sting.utils.io.FileExtension
import java.io.{IOException, PrintWriter, File}
import scala.Enumeration

import org.broadinstitute.sting.utils.exceptions.UserException
import org.apache.commons.io.IOUtils
import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline.Input
import org.broadinstitute.sting.commandline.Output
import scala.io.Source
import org.broadinstitute.sting.utils.exceptions.UserException.CouldNotReadInputFile

import org.broadinstitute.sting.queue.util.StringFileConversions._

import org.broadinstitute.sting.queue.util.Logging

/**
 * Enumeration with different classes of bam files.
 */
object BamType extends Enumeration {
    type BamType = Value
    val TUMOR, NORMAL, SPIKED = Value
}

import BamType._

class AbbreviatedFile( file: File , val abbreviation: String) extends File(file) with FileExtension {
    def withPath(path: String) = new AbbreviatedFile(path, abbreviation)
    
}

class AnnotatedBamFile(file: File , val typeOfBam: BamType, name: String) extends AbbreviatedFile(file, name) {
}

object AnnotatedBamFile extends Logging{
    import BamType._
    def readBamTypesFile(typesFile: File) = {
        try{

            val source = Source.fromFile(typesFile)

            val lines = source.getLines()
            val splits = lines.map(_.split('\t')).map(splitLine => (splitLine(0), splitLine(1), splitLine(2))).toList


            val tumor = splits.collect  { case ("TUMOR", file, name) => new AnnotatedBamFile(file, TUMOR,name)  }
            val normal = splits.collect { case ("NORMAL", file, name) => new AnnotatedBamFile(file, NORMAL,name)  }
            val spiked = splits.collect { case ("SPIKED", file, name) => new AnnotatedBamFile(file, SPIKED,name)  }

            logger.debug("Normal: "+normal)
            logger.debug("Tumor: "+tumor)
            logger.debug("Spiked: " + spiked)

            (normal, tumor, spiked)
        } catch {
            case e: IOException => throw new CouldNotReadInputFile(typesFile, "Could not read bam types file", e)
        }
    }

    def printBamInfoFile(outputFile: File, bamFiles: Seq[AnnotatedBamFile])={
        var writer: PrintWriter = null
        try{
            writer = new PrintWriter(outputFile)
            bamFiles.foreach(file => writer.println( "%s\t%s\t%s".format(file.typeOfBam,file.getAbsolutePath, file.abbreviation) ) )
        } catch {
            case e: IOException => throw new UserException.CouldNotCreateOutputFile(outputFile, "Could not write bam types file.",e)
        } finally {
            IOUtils.closeQuietly(writer)
        }
    }


}


class OutputBamTypeFile extends InProcessFunction{

    @Input
    var bams: Seq[AnnotatedBamFile] = Nil

    @Output
    var bamTypeFile: File = _


    def run() = {
        AnnotatedBamFile.printBamInfoFile(bamTypeFile, bams)
    }

}