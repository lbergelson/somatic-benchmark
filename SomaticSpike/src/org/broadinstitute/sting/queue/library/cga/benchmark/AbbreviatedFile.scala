package org.broadinstitute.sting.queue.library.cga.benchmark

import org.broadinstitute.sting.utils.io.FileExtension
import java.io.{IOException, PrintWriter, File}

import org.broadinstitute.sting.utils.exceptions.UserException
import org.apache.commons.io.IOUtils
import org.broadinstitute.sting.queue.function.InProcessFunction
import org.broadinstitute.sting.commandline.Input
import org.broadinstitute.sting.commandline.Output
import scala.io.Source
import org.broadinstitute.sting.utils.exceptions.UserException.CouldNotReadInputFile

import org.broadinstitute.sting.queue.util.StringFileConversions._
import org.broadinstitute.sting.queue.util.Logging



import BamType._

class AbbreviatedFile( file: File , val abbreviation: String) extends File(file) with FileExtension {
    def withPath(path: String) = new AbbreviatedFile(path, abbreviation)
    
}

class AnnotatedBamFile(file: File , val typeOfBam: BamType, abbreviation: String) extends AbbreviatedFile(file, abbreviation) {
    def toOutputString =  "%s\t%s\t%s".format(typeOfBam,getAbsolutePath, abbreviation)
}

object AnnotatedBamFile extends Logging{
    class AnnotatedBamException extends Exception

    def apply(str: String) = {
        val splitStr = str.split('\t')
        try {
            new AnnotatedBamFile(splitStr(0), BamType.withName(splitStr(1)), splitStr(2) )
        } catch {
            case e: Exception => throw new AnnotatedBamException()
        }
    }

    def readBamTypesFile(typesFile: File) = {
        try{
            val source = Source.fromFile(typesFile)
            val lines = source.getLines()
            val annotated = lines.map(AnnotatedBamFile(_)).toSeq

            warnIfAbbreviationsAreNotUnique(annotated)

            annotated
        } catch {
            case e: IOException => throw new CouldNotReadInputFile(typesFile, "Could not read bam types file due to IO problem:", e)
            case e: AnnotatedBamException => throw new CouldNotReadInputFile(typesFile, "Couldn't load the bam types file because of a problem creating an AnnotatedBamFile", e)
        }
    }

    def splitByType(files: Seq[AnnotatedBamFile]) = {
        val tumor = files.filter( _.typeOfBam == TUMOR )
        val normal = files.filter( _.typeOfBam == NORMAL )
        val spiked = files.filter(_.typeOfBam == SPIKED)

        logger.debug("Normal: "+normal)
        logger.debug("Tumor: "+tumor)
        logger.debug("Spiked: " + spiked)

        (normal, tumor, spiked)
    }

    def printBamInfoFile(outputFile: File, bamFiles: Seq[AnnotatedBamFile])={
        warnIfAbbreviationsAreNotUnique(bamFiles)
        var writer: PrintWriter = null
        try{
            writer = new PrintWriter(outputFile)
            bamFiles.foreach(file => writer.println(file.toOutputString ) )
        } catch {
            case e: IOException => throw new UserException.CouldNotCreateOutputFile(outputFile, "Could not write bam types file.",e)
        } finally {
            IOUtils.closeQuietly(writer)
        }
    }

    private def namesAreUnique(bamFiles: Seq[AnnotatedBamFile]) = {
        bamFiles.map(_.abbreviation).distinct.length != bamFiles.length
    }

    private def warnIfAbbreviationsAreNotUnique(bamFiles: Seq[AnnotatedBamFile]) = {
        if(!namesAreUnique(bamFiles)) {
            logger.warn("Bam file abbreviations are not unique!  This will probably cause problems!")
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