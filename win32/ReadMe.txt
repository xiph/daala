/*Daala video codec
Copyright (c) 2001-2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

Before you can build this Daala source code on a Windows machine,
you have to obtain and build a few required dependencies.

=== Dependencies ===
* http://www.visualstudio.com/en-us/downloads/download-visual-studio-vs.aspx
  Visual Studio versions 2008 or 2010 are recommended
* http://downloads.xiph.org/releases/ogg/?C=M;O=D
  libogg v1.3.0 or later is recommended. The VS solutions have been tested
  using libogg 1.3.2: http://downloads.xiph.org/releases/ogg/libogg-1.3.2.zip
* https://www.wxwidgets.org/downloads/
  wxWidgets is only required if you want to run the Daala bitstream analyser.
  The VS solutions have been tested with version 3.0.2 from
  http://sourceforge.net/projects/wxwindows/files/3.0.2/)

The VS solution files were created using Visual Studio 2008 Team System
and Visual Studio 2010 Ultimate. Even though the provided solutions might work
with the express versions of VS, this has not been tested.

=== Installation ===
Clone the Daala.git repository into a folder on your machine (this can take
several minutes). The clone operation should by default create a new "daala"
folder.

    git clone https://git.xiph.org/daala.git

Make sure you run the git clone operation on the same machine where you intend
to use the code. Checking out a copy on Windows and then trying to use it on
Linux will not work, as executable permissions and line-endings will not be
set properly.

Unpack the previously downloaded libogg into a folder named "ogg" in the same
folder you cloned daala into.

If you also want to be able to build the Stream Analyser, you need to:
- make a "wxWidgets" folder in the same folder your "daala" folder is
- either Git clone/fork or download wxWidgets into it
- within wxWidgets\build\msw, open the relevant MSVC solution and build the
  Debug and/or Release versions of wxWidgets, depending on if you need to run
  the Daala executables in debug-mode or normally.

Your folder structure should now look like this:

    <your folder>\daala
    <your folder>\ogg
    <your folder>\wxWidgets

== Building the executables ==

* Open the file "Daala_static.sln" located within the relevant
  "daala\win32\Visual Studio\VS20xx" folder.
* Build the solution.

This will generate 3 executables and 3 static libraries:

* encoder_example.exe
* encoder_example.exe
* decoder_example.exe
* LibDaalaBase.lib
* libdaaladec.lib
* libdaalaenc.lib

You can find more information at https://wiki.xiph.org/Daala or on the #daala
channel on irc.freenode.net.