/*
  FILENAME:   SampleECG.ino
  AUTHOR:     Orlando S. Hoilett
  EMAIL:      orlandohoilett@gmail.com
  VERSION:    0.0

  
  AFFILIATIONS:
  BME 695: Measurements Fall 2017
  (Purdue University, West Lafayette, IN)


  UPDATES:
  Version 0.0
  2017/12/17:1755>
              First version of the code.

              
  DESCRIPTION  
  This code uses the Arduino as sipmly a serial interface between
  my ECG analog front-end and MATLAB. This code was developed as
  part of an extra credit project for BME 695: Measurements
  (Fall 2017) class.

  The code uses timer interrupts to achieve a very consistent
  sampling rate for sampling the ECG signal.

    
  DISCLAIMER
  This code is in the public domain. Please feel free to modify,
  use, etc however you see fit. But, please give reference to
  original authors as a courtesy to Open Source developers.

*/

const uint8_t ecg = A0; //analog input pin

void setup ()
{
  Serial.begin(115200);
  analogReference(EXTERNAL);
  
  //set timer1 interrupt at 1kHz
  cli(); //first disable interrupts in order to configure them
  TCCR1A = 0; // set entire TCCR1A register to 0
  TCCR1B = 0; // same for TCCR1B
  TCNT1  = 0; //initialize counter value to 0

  //compare match register =
  //[ MCLK (Hz)/ (prescaler * fsample (Hz)) ] - 1
  //MCLK for Arduino Uno is 16MHz
  //fsample is desired sampling frequency in Hertz
  //use 15 for Fs = 1041Hz
  //use 78 for Fs = 200Hz
  //use 156 for Fs = 100Hz
  //use 260 for Fs = 60Hz
  //compare match register = [ 16MHz/ (prescaler * desired interrupt frequency) ] - 1
  OCR1A = 78;
  
  TCCR1B |= (1 << WGM12); // turn on CTC mode
  TCCR1B |= (1 << CS12) | (1 << CS10); // Set CS10 and CS12 bits for 1024 prescaler
  TIMSK1 |= (1 << OCIE1A); // enable timer compare interrupt
  sei(); //now that interrupts have been set, enable them
}


void loop()
{
  //everything is handled in the interrupt service routine (ISR)
}


//Interrupt Service Routine
//Samples the ECG signal from a selectd analog input and prints
//the data to the serial port. The frequency of the ISR is set
//in the setup() function
ISR(TIMER1_COMPA_vect)
{
  Serial.println(analogRead(ecg));
}
