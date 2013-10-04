package org.broadinstitute.sting.queue.library.cga.benchmark

/**
 * Enumeration with different classes of bam files.
 */
object BamType extends Enumeration {
    type BamType = Value
    val TUMOR, NORMAL, SPIKED = Value
}
