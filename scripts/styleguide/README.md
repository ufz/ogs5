OpenGeoSys C++ Style Guide
==========================

This styleguide is derived from the [Google C++ Style Guide][googlestyle]. Some modifications were made to match the [OpenGeoSys][ogs] development needs.

Automatically apply this guide
------------------------------

This guide comes with a [configuration file][config] for the [uncrustify][crust] automatic code styler. With this tool you can restyle your source code according to this styleguide. The config file is not finished yet. For creating the config file the [UniversalIndentGUI][gui] was used.

To reformat your code:

- Install [uncrustify][crust]
- Download the [configuration file][config]
- Download the [uncrustify.sh script][script]
- Run the script on your code:

		uncrustify.sh [Path to your code directory] [file extension]
		E.g for header files: uncrustify.sh sources/ h
		E.g for source files: uncrustify.sh sources/ cpp


[googlestyle]: http://code.google.com/p/google-styleguide/
[ogs]: http://www.opengeosys.net/
[config]: http://ufz.github.com/styleguide/uncrustify.cfg
[crust]: http://uncrustify.sourceforge.net/
[gui]: http://universalindent.sourceforge.net/
[script]: http://ufz.github.com/styleguide/uncrustify.sh