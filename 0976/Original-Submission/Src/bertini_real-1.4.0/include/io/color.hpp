
#ifndef BERTINI_REAL_COLOR_HPP
#define BERTINI_REAL_COLOR_HPP

#include <string>




/**
 \brief terminal-control of colors.
 
 namespace of color functions, using terminal controls.
 
 e.g. \033[0;30m for black.
 */
namespace color {
	
	/** get the integer associated with the name of a color 
	 
	 k - black - 30
	 r - red - 31
	 g - green - 32
	 y - brown - 33
	 b - blue - 34
	 m - magenta - 35
	 c - cyan - 36
	 l - lightgray - 37
	 
	 \return integer corresponding to color
	 \param c the single-character string of the color.
	 */
	std::string color_to_int(const char c);
	
	/**
	 set the text to bold
	 
	 \033[1;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string bold(char new_color);
	
	/**
	 set the text to darker color
	 
	 \033[2;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string dark(char new_color);
	
	/**
	 set the text to underline
	 
	 \033[4;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string underline(char new_color);
	
	/**
	 set the background to a color indicated by a name
	 
	 \033[7;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string background(char new_color);
	
	/**
	 set the text to strikethrough
	 
	 \033[9;XXm
	 
	 \return string to print to cout to control color
	 \param new_color the single-character name of the color
	 */
	std::string strike(char new_color);
			
	
	
	
	/**
	 set the text to whatever the console believes is the default
	 
	 \033[0m
	 
	 \return string to print to cout to control color
	 */
	std::string console_default();
	
	/**
	 set the text to black
	 
	 \033[0;30m
	 
	 \return string to print to cout to control color
	 */
	std::string black();
	
	/**
	 set the text to red
	 
	 \033[0;31m
	 
	 \return string to print to cout to control color
	 */
	std::string red();
	
	/**
	 set the text to green
	 
	 \033[0;32m
	 
	 \return string to print to cout to control color
	 */
	std::string green();

	
	/**
	 set the text to brown
	 
	 \033[0;33m
	 
	 \return string to print to cout to control color
	 */
	std::string brown();
	
	/**
	 set the text to blue
	 
	 \033[0;34m
	 
	 \return string to print to cout to control color
	 */
	std::string blue();
	
	/**
	 set the text to magenta
	 
	 \033[0;35m
	 
	 \return string to print to cout to control color
	 */
	std::string magenta();
	
	/**
	 set the text to cyan
	 
	 \033[0;36m
	 
	 \return string to print to cout to control color
	 */
	std::string cyan();
	
	/**
	 set the text to gray
	 
	 \033[0;37m
	 
	 \return string to print to cout to control color
	 */
	std::string gray();
	
	
	//black - 30
	//red - 31
	//green - 32
	//brown - 33
	//blue - 34
	//magenta - 35
	//cyan - 36
	//lightgray - 37
	
}


#endif

