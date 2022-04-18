TEMPLATE      = app
CONFIG        = warn_on release exceptions

TARGET        = precompile

DEFINES       = #_DEBUG_

INCLUDEPATH   = include \
                include/xml 

OBJECTS_DIR   = obj

MOC_DIR       = $$OBJECTS_DIR

HEADERS      += include/xml/attribute.h \
                include/xml/commentnode.h \
                include/xml/elementnode.h \
                include/xml/pinode.h \
                include/xml/textnode.h \
                include/xml/xmlnode.h \
                include/xml/XMLParseFunctions.h \
                include/error.h \
                include/ConfigElem.h \
                include/ConfigFile.h \
                include/Convert.h \
                include/defs.h \
                include/ParseCmdLine.h \
                include/ParseElem.h \
                include/ParseQueue.h \
                include/Parser.h \
                include/ParserConf.h \
                include/PreParser.h \
                include/PrintArg.h \
                include/PrintConverter.h \
                include/Reserved.h \
                include/Var.h \
                include/Variables.h \
                include/ConvertElem.h \
                include/MatchElem.h \
                include/Matcher.h \
                include/Graph.h \
                include/ConvertConfig.h \
                include/finder.h \
                include/xml/XMLParser.h 

SOURCES      += src/xml/attribute.cpp \
                src/xml/commentnode.cpp \
                src/xml/elementnode.cpp \
                src/xml/pinode.cpp \
                src/xml/textnode.cpp \
                src/xml/xmlnode.cpp \
                src/ConfigElem.cpp \
                src/ConfigFile.cpp \
                src/Convert.cpp \
                src/ParseCmdLine.cpp \
                src/ParseElem.cpp \
                src/ParseQueue.cpp \
                src/Parser.cpp \
                src/ParserConf.cpp \
                src/precompile.cpp \
                src/PreParser.cpp \
                src/PrintArg.cpp \
                src/PrintConverter.cpp \
                src/Reserved.cpp \
                src/Var.cpp \
                src/Variables.cpp \
                src/ConvertElem.cpp \
                src/MatchElem.cpp \
                src/Matcher.cpp \
                src/Graph.cpp \
                src/ConvertConfig.cpp \
                src/finder.cpp \
                src/xml/XMLLexer.cpp \
                src/xml/XMLParser.cpp \
                src/lex.yy.cpp
