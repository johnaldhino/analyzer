###### Hall A Software hana_decode SConscript Build File #####
###### Author:  Edward Brash (brash@jlab.org) June 2013

import os
import re
import SCons.Util
Import('baseenv')

standalone = baseenv.subst('$STANDALONE')
print ('Compiling decoder executables:  STANDALONE = %s\n' % standalone)

standalonelist = Split("""
tstoo tstfadc tstf1tdc tstio tdecpr prfact epicsd tdecex tst1190
""")
# Still to come, perhaps, are (etclient, tstcoda) which should be compiled
# if the ONLINE_ET variable is set.  

list = Split("""
Caen1190Module.C
Caen775Module.C
Caen792Module.C
CodaDecoder.C
F1TDCModule.C
Fadc250Module.C
FastbusModule.C
GenScaler.C
Lecroy1875Module.C
Lecroy1877Module.C
Lecroy1881Module.C
Module.C
PipeliningModule.C
Scaler1151.C
Scaler3800.C
Scaler3801.C
Scaler560.C
SimDecoder.C
THaCodaData.C
THaCodaDecoder.C
THaCodaFile.C
THaCrateMap.C
THaEpics.C
THaEvData.C
THaFastBusWord.C
THaSlotData.C
THaUsrstrutils.C
VmeModule.C
""")

# Requires SCons >= 2.3.5 for "exclude" keyword
#list = Glob('*.C',exclude=['*_main.C','*_onl.C','calc_thresh.C',
#                           'THaEtClient.C','THaGenDetTest.C'])

#baseenv.Append(LIBPATH=['$HA_DIR'])

proceed = "1" or "y" or "yes" or "Yes" or "Y"
if baseenv.GetOption('clean'):
	print ('!!!!!! Cleaning hana_decode standalone executables !!!!!!')
	for scalex in standalonelist:
		hana_decode_clean_command = "rm %s; rm ../bin/%s" % (scalex,scalex)
		print ("hana_decode_clean_command = %s" % hana_decode_clean_command)
		os.system(hana_decode_clean_command)

if baseenv.subst('$STANDALONE')==proceed:
        for scalex in standalonelist:
                pname = scalex

		if scalex=='epicsd':
			main = 'epics_main.C'
		else:
                	main = scalex+'_main.C'

                if scalex=='tdecex':
			pname = baseenv.Program(target = pname, source = [main,'THaGenDetTest.C'])
		else:	
			pname = baseenv.Program(target = pname, source = [main])
			
                baseenv.Install('../bin',pname)
                baseenv.Alias('install',['../bin'])

sotarget = 'dc'

#dclib = baseenv.SharedLibrary(target=sotarget, source = list+['THaDecDict.so'],SHLIBPREFIX='../lib',SHLIBVERSION=['$VERSION'],LIBS=[''])
dclib = baseenv.SharedLibrary(target=sotarget, source = list+['THaDecDict.so'],SHLIBPREFIX='../lib',LIBS=[''])
print ('Decoder shared library = %s\n' % dclib)

linkbase = baseenv.subst('$SHLIBPREFIX')+sotarget

cleantarget = linkbase+'.so.'+baseenv.subst('$VERSION')
majorcleantarget = linkbase+'.so'
localmajorcleantarget = '../'+linkbase+'.so'
shortcleantarget = linkbase+'.so.'+baseenv.subst('$SOVERSION')
localshortcleantarget = '../'+linkbase+'.so.'+baseenv.subst('$SOVERSION')

print ('cleantarget = %s\n' % cleantarget)
print ('majorcleantarget = %s\n' % majorcleantarget)
print ('shortcleantarget = %s\n' % shortcleantarget)
try:
	os.symlink(cleantarget,localshortcleantarget)
	os.symlink(shortcleantarget,localmajorcleantarget)
except:	
	print (" Continuing ... ")

Clean(dclib,cleantarget)
Clean(dclib,localmajorcleantarget)
Clean(dclib,localshortcleantarget)

#baseenv.Install('../lib',dclib)
#baseenv.Alias('install',['../lib'])
