The reository acts in two way;

 1. Everything which is named as tutorial and/or Intro are intended for the people who want to learn radio astronomy, it mostly covers the explanation from ADC data to interferometry and detection of sun, galactic plane, and all type of solar burst in one tutorial code. Pulsar tutorial will consist of understaiding the data types used at different level of data proccessing like - .PACAP, .DADA, .FIL, .FITS and uses diffferent python libraries to understand them along with introduction of psrchive;


### Every codes are written from scrath in tutorial codes, so it also helps sometimes if you want to play with the data for understadig; which i keep on doing.

 2. It also has the 'def' (python refrence) and source codes for all the data proccesing pipelines which i have written for the following telescopes :

     1. two-element interferometry which it has been used and currently using for the study of solar and            gaalctic transits along with the dailyobservations of solar flares (detection of solar flares was not done by me.)
     2. Diamond Array Telescope, currently writing the transient detection pipeline which consist of the whole data proccesing from start of ADC data products to classifications of the same; parralell proccessing for ADC to filterbank by modifying the old codes, and gpu and/or search algorithms for transients along with rfi mitigation. (The codes order may look messy but it will be better as the time goes hopefully by the time we it a real time transient detection system.)
