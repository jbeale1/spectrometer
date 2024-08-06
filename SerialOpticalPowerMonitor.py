#!/usr/bin/python3
# Python code to read data via USB serial port on "FieldBest" optical power meter
# 2024-July-16 J.Beale

from serial import *
import time, datetime

LOGFILE = "PowerLog.csv"
PORT = "/dev/ttyACM0"
VERSION  = "Optical Power Meter 7/16/2024 JPB"
CSV_HEADER = "date_time, epoch, power"
eol_str = "\n"  # end of line string in file output

timeInterval = 1   # how many seconds between readings

# at a rate of 5 samples per second, every 4th reading is a duplicate
# 3 dups in 10 samples (2 seconds)
ser=Serial(PORT,19200,8,'N',1,timeout=0.3)  # serial input

dstring = datetime.datetime.today().strftime('%Y-%m-%d_%H%M%S_')

f = open(dstring+LOGFILE, 'w')  # open log file
print("%s" % VERSION)
print("%s" % CSV_HEADER)
f.write ("%s\n" % CSV_HEADER)
f.write ("# %s\n" % VERSION)

outcmd = '?pw%\n'.encode('utf_8')  # encode power meter command string as bytes

ser.read(255)            # clear out existing buffer & throw it away
lines = 0                # how many lines we've read
while True:
    buf = ""
    t = datetime.datetime.now()
    """
    sleeptime = timeInterval - ((t.second % timeInterval) + t.microsecond/1000000.0)
    tStop = t + datetime.timedelta(seconds=sleeptime)
    while True:      # wait until the next even 10 second mark
      ser.write(outcmd)
      buf = ser.readline()              # byte string
      tEpoch = int(time.time())         # time() = seconds.usec since Unix epoch
      tNow = datetime.datetime.now()
      if (tNow > tStop):
          break
    """
    ser.write(outcmd)       # write this byte string command to device
    buf = ser.readline()    # get response

    tEpoch = (int(time.time()*10))/10.0         # time() = seconds.usec since Unix epoch
    tNow = datetime.datetime.now()

    # print(buf)
    ptr = buf.find(b'\x00') # 0x00 is used as a separator character
    power = buf[0:ptr].decode("utf-8") # optical power in mW
    scale = buf[ptr+1:].decode("utf-8") # some sort of scaling factor
    # buffer = ("%s,%s" % (power, scale))
    buffer = ("%s" % (power))
    
    # buffer = buf.decode("utf-8").rstrip()         # string terminated with '\n'
    if (buffer != '') :
      outbuf = str(datetime.datetime.now())[:-3] + "," + \
        ("%12.1f" % tEpoch ) + "," + buffer
      print (outbuf)
      f.write (outbuf)
      f.write (eol_str)
      lines += 1
      if (lines % 10 == 0):  # save it out to actual file, every so often
          f.flush()

f.close                  # close log file
ser.close()            # close serial port when done. If we ever are...
