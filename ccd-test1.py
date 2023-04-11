#!/usr/bin/env python
""" 
communicates with this Toshiba TCD1304 ccd + STM32F103C8T device, a black PCB marked "YX_TCD1304 ver1.7"
from www.aliexpress.us/item/3256803138125304.html
On linux, needs https://www.wch.cn/downloads/CH341SER_LINUX_ZIP.html
because as of April 2023, mainline linux CH341 driver is old & broken at this bitrate
  J.Beale, April 10 2023
"""

"""
0xA1   Request 12-bit data command
Return value: 7296 bytes. 3648 pixels, the lower eight-bit data is in front, and the upper eight-bit data is in the back)

0xB1     integration time 10 us
0xB2     integration time 20 us
0xB3     integration time 50 us
0xB4     integration time 60 us
0xB5     integration time 75 us
0xB6     integration time 100 us
0xB7     integration time 500 us
0xB8     integration time 1.25 ms
0xB9     integration time 2.5 ms  
0xBA     integration time 7.5 ms   
"""

import serial
import time
import sys

serport = "/dev/ttyCH341USB0"
fnameout = "/home/john/ccd-dat1.pgm" # PGM = portable greymap

xsize = 3648 # pixels in linear CCD array
ysize = 500  # how many lines to read

def readSerPort(ser):
	while (ser.in_waiting > 0):
		bytesIn = ser.in_waiting
		serdata = ser.read(bytesIn)
		print(serdata.hex())


# -----------------------------------------
ser = serial.Serial(serport, baudrate=921600)  # module uses this baud rate

f = open(fnameout, 'wb')
header = ("P5 %d %d 255\n" % (xsize, ysize)).encode('utf_8')
f.write(header)

# --------------------------------------------------
# packet = bytes([0xBA])  # 7.5 msec (longest integration time)
packet = bytes([0xB3])  # set integration time
ser.write(packet)
time.sleep(0.1)  # need some delay, but not sure how much

packet = bytearray()
#packet.append(0xA1)  # read TCD1304 module (12 bpp, 7296 bytes)
packet.append(0xA2)  # read TCD1304 module (8 bpp, 3648 bytes)
# packet.append(0xA3)  # read data min/max summary from firmware

for i in range(ysize):
	ser.write(packet)
	time.sleep(0.08)
	rawdata = ser.read(ser.in_waiting)
	bcount = len(rawdata)
	if (bcount != xsize):
		print("Error: pixel count = %d" % bcount)
		sys.exit()
	# print(bcount)
	f.write(rawdata)


ser.close()  # close the port
f.close()    # close the output file
