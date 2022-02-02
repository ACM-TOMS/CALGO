#include "io/color.hpp"

namespace color {
	std::string color_to_int(const char c)
	{
		switch (c) {
			case 'k':
				return "30";
			case 'r':
				return "31";
			case 'g':
				return "32";
			case 'y':
				return "33";
			case 'b':
				return "34";
			case 'm':
				return "35";
			case 'c':
				return "36";
			case 'l':
				return "37";
			default:
				return "30";
		}
	}
	
	
	std::string bold(char new_color)
	{
		return "\033[1;" + color_to_int(new_color) + "m";
	}
	
	std::string dark(char new_color)
	{
		return "\033[2;" + color_to_int(new_color) + "m";
	}
	
	
	std::string underline(char new_color)
	{
		return "\033[4;" + color_to_int(new_color) + "m";
	}
	
	
	std::string background(char new_color)
	{
		return "\033[7;" + color_to_int(new_color) + "m";
	}
	
	
	std::string strike(char new_color)
	{
		return "\033[9;" + color_to_int(new_color) + "m";
	}
	
	std::string console_default(){
		return "\033[0m";
	}
	
	std::string black(){
		return "\033[0;30m";
	}
	
	std::string red(){
		return "\033[0;31m";
	}
	
	std::string green(){
		return "\033[0;32m";
	}
	
	
	std::string brown(){
		return "\033[0;33m";
	}
	
	std::string blue(){
		return "\033[0;34m";
	}
	
	std::string magenta(){
		return "\033[0;35m";
	}
	
	std::string cyan(){
		return "\033[0;36m";
	}
	
	std::string gray(){
		return "\033[0;37m";
	}
	
	
	//black - 30
	//red - 31
	//green - 32
	//brown - 33
	//blue - 34
	//magenta - 35
	//cyan - 36
	//lightgray - 37
}

