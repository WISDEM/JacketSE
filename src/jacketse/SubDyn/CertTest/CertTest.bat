@ECHO OFF
@ECHO.


REM  Set up environment variables.  You will probably have to change these.


@SET Compare=FC /T
@SET CRUNCH=..\bin\crunch_win32.exe
rem @SET CRUNCH=Call Crunch

@SET MATLAB=matlab
@SET DateTime=DateTime.exe
@SET Editor=NotePad.EXE
@SET CompareFile=CertTest.out

@SET ResultsDir=NREL_RESULTS

::=======================================================================================================
IF /I "%1"=="-DEBUG" GOTO debugVer
IF /I "%1"=="-GFORTRAN" GOTO gfortran
IF /I "%1"=="-IFORT" GOTO ifort
IF /I "%1"=="-DEVBUILD" GOTO devBuild
IF /I "%1"=="-DEVDEBUG" GOTO devDebugBuild

:releaseVer
@SET EXE_VER=Using released version of SubDyn (IVF/VS)
@SET SubDyn=..\bin\SubDyn_win32.exe
goto CertTest

:debugVer
@SET EXE_VER=Using SubDyn compiled in debug mode (IVF/VS)
@SET SubDyn=..\bin\SubDyn_debug_win32.exe
goto CertTest

:gfortran
@SET EXE_VER=Using SubDyn compiled with makefile (gfortran)
@SET SubDyn=..\compiling\SubDyn_gwin32.exe
goto CertTest

:ifort
@SET EXE_VER=Using SubDyn compiled with Compile_SubDyn.bat (IVF)
@SET SubDyn=..\compiling\SubDyn_iwin32.exe
goto CertTest

:devBuild
@SET EXE_VER=Using SubDyn compiled with Visual Studio Project, release mode (IVF/VS)
@SET SubDyn=..\bin\SubDyn_dev_Release_win32.exe
goto CertTest

:devDebugBuild
@SET EXE_VER=Using SubDyn compiled with Visual Studio Project, debug mode (IVF/VS)
@SET SubDyn=..\bin\SubDyn_dev_debug_win32.exe
goto CertTest

::=======================================================================================================


:CertTest


REM  SubDyn test sequence definition:

@SET  TEST01=Test #01: Version of the OC3 Monopile.
@SET  TEST02=Test #02: A 2D house frame aligned along the X-axis.
@SET  TEST03=Test #03: A one-bay OC4 Jacket under static displacement.
@SET  TEST04=Test #04: OC4 Jacket under static displacement.
@SET  TEST05=Test #05: A 2D house frame rotated 30 deg from the X-axis.

@SET  DASHES=---------------------------------------------------------------------------------------------
@SET  POUNDS=#############################################################################################

@IF EXIST CertTest.out  DEL CertTest.out

ECHO.                                               >> CertTest.out
ECHO           ************************************ >> CertTest.out
ECHO           ** SubDyn Acceptance Test Results ** >> CertTest.out
ECHO           ************************************ >> CertTest.out

ECHO.                                                                             >> CertTest.out
ECHO ############################################################################ >> CertTest.out
ECHO # Inspect this file for any differences between your results and the saved # >> CertTest.out
ECHO # results.  Any differing lines and the two lines surrounding them will be # >> CertTest.out
ECHO # listed.  The only differences should be the time stamps at the start of  # >> CertTest.out
ECHO # each file.                                                               # >> CertTest.out
ECHO #                                                                          # >> CertTest.out
ECHO # If you are running on something other than a PC, you may see differences # >> CertTest.out
ECHO # in the last significant digit of many of the numbers.                    # >> CertTest.out
ECHO ############################################################################ >> CertTest.out

rem ECHO.                                            >> CertTest.out
rem ECHO Date and time this acceptance test was run: >> CertTest.out
rem %DateTime%                                       >> CertTest.out
rem ECHO.                                            >> CertTest.out


ECHO.                                            >> CertTest.out
ECHO %EXE_VER%                                   >> CertTest.out
ECHO SubDyn called with this command:              >> CertTest.out
ECHO %SubDyn%                                      >> CertTest.out
ECHO.                                            >> CertTest.out


echo %DASHES%
echo %EXE_VER%
echo %SubDyn%
echo %DASHES%


rem *******************************************************
:Test1
@CALL :GenTestHeader %Test01%
@CALL :RunSubDyn 01 out

rem *******************************************************
:Test2
@CALL :GenTestHeader %Test02%
@CALL :RunSubDyn 02 out

rem *******************************************************
:Test3
@CALL :GenTestHeader %Test03%
@CALL :RunSubDyn 03 out

rem *******************************************************
:Test4
@CALL :GenTestHeader %Test04%
@CALL :RunSubDyn 04 out

rem *******************************************************
:Test5
@CALL :GenTestHeader %Test05%
@CALL :RunSubDyn 05 out

rem ******************************************************
rem  Let's look at the comparisons.
:MatlabComparisons

rem %MATLAB% /r PlotCertTestResults('.','.\TstFiles');exit;


rem %Editor% CertTest.out
goto END

rem ******************************************************
:GenTestHeader
echo %POUNDS%
@echo SubDyn %*
echo %POUNDS%

echo.                                    >> %CompareFile%
echo %POUNDS%                            >> %CompareFile%
echo.                                    >> %CompareFile%
echo %*                                  >> %CompareFile%
EXIT /B

rem ******************************************************
:RunSubDyn
:: Run SubDyn.
@SET TEST=%1

%SubDyn% Test%1.dvr


IF ERRORLEVEL 1  GOTO ERROR

rem echo %DASHES%                          >> %CompareFile%
rem %Compare% Test%1\Test%1.SD.out Test%1\%ResultsDir%\Test%1.SD.out >> %CompareFile%
rem %Compare% Test%1\Test%1.SD.sum Test%1\%ResultsDir%\Test%1.SD.sum >> %CompareFile%


echo %DASHES%



EXIT /B




:ERROR
:: Sets clears memory and stops the batch immediately
@echo ** An error has occurred in Test #%TEST% **
@echo ** An error has occurred in Test #%TEST% ** >> %CompareFile%

@call :end
@call :__ErrorExit 2> nul
EXIT /B

:__ErrorExit
rem Creates a syntax error, stops immediately
()
EXIT /B


:END

@SET CRUNCH=
@SET MATLAB=
@SET MBC_SOURCE=
@SET Compare=
@SET CompareFile=
@SET DASHES=
@SET DateTime=
@SET Editor=
@SET SubDyn=
@SET POUNDS=
@SET TEST=
@SET TEST01=
@SET TEST02=
@SET TEST03=
@SET TEST04=


SET EXE_VER=

rem type Bell.txt
@echo Processing complete.

EXIT /B