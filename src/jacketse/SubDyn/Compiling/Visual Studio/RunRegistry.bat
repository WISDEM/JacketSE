@ECHO OFF

set lines=======================================================================
echo %lines%
IF "%1"=="" (
ECHO.
ECHO   The calling syntax for this script is
ECHO             RunRegistry ModuleName
ECHO.
GOTO Done
)


REM ----------------------------------------------------------------------------
REM ------------------------- LOCAL PATHS --------------------------------------
REM ----------------------------------------------------------------------------
REM -- USERS MAY EDIT THESE PATHS TO POINT TO FOLDERS ON THEIR LOCAL MACHINES. -
REM -- NOTE: do not use quotation marks around the path names!!!! --------------
REM ----------------------------------------------------------------------------
REM ----------------------------------------------------------------------------

SET Registry=..\..\bin\Registry_win32.exe
SET Source_Loc=..\..\Source

SET NWTC_Lib_Loc=%Source_Loc%\dependencies\NWTC_Library





IF /I "%2"=="dev" CALL ..\Set_FAST_paths.bat



REM ----------------------------------------------------------------------------
REM ---------------- RUN THE REGISTRY TO AUTO-GENERATE FILES -------------------
REM ----------------------------------------------------------------------------


SET CURR_LOC=%Source_Loc%
%REGISTRY% "%CURR_LOC%\SubDyn_Registry.txt" -I %NWTC_Lib_Loc%
GOTO checkError



:checkError
ECHO.
IF %ERRORLEVEL% NEQ 0 (
ECHO Error running  Registry for SubDyn.
) ELSE (
ECHO SubDyn_Types.f90 was created.
COPY /Y "SubDyn_Types.f90" "%CURR_LOC%"
)




:end
REM ----------------------------------------------------------------------------
REM ------------------------- CLEAR MEMORY -------------------------------------
REM ----------------------------------------------------------------------------
ECHO. 


SET REGISTRY=

SET NWTC_Lib_Loc=
SET Source_Loc=

SET CURR_LOC=
:Done
echo %lines%
set lines=