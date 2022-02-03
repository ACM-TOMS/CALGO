package com.googlecode.compressingcircularbuffers;
/** Compresses by aggregating adjacent samples into their totals (sum).*/
public class SummingCompressionMethod implements CompressionMethod {
  private double sum = 0;
  private int count = 0;
  public void update(double value) {
    sum += value;
    count++;
  }
  public void reset() {
    sum = 0;
    count = 0;
  }
  public double getCompressedValue() {
    return sum;
  }
  public double mergeValues(double older, double newer) {
    double result = older + newer;
    return result;
  }
  public int getCount() { return count; }
}
