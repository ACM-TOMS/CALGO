#pragma once

#include <vector>
/**
 \brief An odometer of odometers, so to speak.
 
 */
class DoubleOdometer
{
	
private:
	
	/**
	 roll one mile.
	 \return 1 if rolled over, 0 if not.  if 1, then need to increment the functions.
	 */
	int increment_registers(){
		
		int carry = 1; // seed carry so that it causes addition of at least the last entry of the odometer
		for (int ii=num_active_registers-1; ii>=0; ii--) { // count down from the end of the indexes
			
			if (carry==1)
				registers[ii]++;
			
			if ( registers[ii]>=(bases[active_registers[ii]]) ) {
				registers[ii] = 0;
				carry = 1;
			}
			else{
				carry = 0;
				break;
			}
		}
		return carry;  // if return 1, then need to increment the functions.
		
	};
	
	
	int increment_active_registers(){
		int carry = 1; // seed 1
		
		for (int ii=num_active_registers-1; ii>=0; ii--) {
			
			if (carry==1){
				
				active_registers[ii]++;
				for (int jj=ii+1; jj<num_active_registers; jj++) {
					active_registers[jj] = active_registers[jj-1]+1;
				}
			}
			
			if (active_registers[num_active_registers-1]>=num_total_registers) {
				carry = 1;
			}
			else{
				carry = 0;
				break; // all done!
			}
			
		}
		
		int local_counter = 0;
		for (int ii=0; (ii<num_total_registers) && (local_counter<num_inactive_registers) ; ii++) {
			if (std::find(active_registers.begin(), active_registers.end(), ii)==active_registers.end())
			{
				inactive_registers[local_counter] = ii;
				local_counter++;
			}
		}
		
		return carry;
	};
	
public:
	
	int num_total_registers;
	int num_active_registers;
	int num_inactive_registers;
	// create and seed the function indices -- keep track of which functions we are working on
	
	std::vector< int > inactive_registers; // of length total - active
	std::vector< int > active_registers; // of length num_active_registers
	
	std::vector< int > bases; // of length num_total_registers
	std::vector< int > registers;
	
	
	/**
	 constructor
	 */
	DoubleOdometer()
	{
		num_total_registers = num_active_registers = num_inactive_registers = 0;
	}
	
	
	/**
	 constructor
	 
	 \param num_total_ the number of total registers there will be.
	 \param num_active_ the number of active registers there will be total.  this must be ≤ num_total_.  If == num_total_, this is a traditional car odometer.
	 \param uniform_base the base of all registers.  if 10, this is a traditional car odometer
	 */
	DoubleOdometer(int num_total_, int num_active_, int uniform_base)
	{
		num_total_registers = num_total_;
		num_active_registers = num_active_;
		num_inactive_registers = num_total_registers - num_active_registers;
		
		for (int ii=0; ii<num_active_registers; ii++)
			active_registers.push_back(ii);
		
		for (int ii=num_active_registers; ii<num_total_registers; ii++)
			inactive_registers.push_back(ii);
		
		for (int ii=0; ii<num_active_registers; ii++)
			registers.push_back(0);
		
		for (int ii=0; ii<num_total_registers; ii++)
			bases.push_back(uniform_base);
		
	}
	
	/**
	 constructor
	 
	 \param num_total_ the number of total registers there will be.
	 \param num_active_ the number of active registers there will be total.  this must be ≤ num_total_.  If == num_total_, this is a traditional car odometer.
	 \param new_bases the base of all registers.  if all 10, this is a traditional car odometer.  The number in here must match num_total_.
	 */
	DoubleOdometer(int num_total_, int num_active_, const std::vector<int> & new_bases)
	{
		num_total_registers = num_total_;
		num_active_registers = num_active_;
		num_inactive_registers = num_total_registers - num_active_registers;
		
		if ( num_total_!= int(new_bases.size()) ) {
			throw std::logic_error("size mismatch in creation of double odometer.  num_total must equal size of base");
		}
		
		for (int ii=0; ii<num_active_registers; ii++)
			active_registers.push_back(ii);
		
		for (int ii=num_active_registers; ii<num_total_registers; ii++)
			inactive_registers.push_back(ii);
		
		for (int ii=0; ii<num_active_registers; ii++)
			registers.push_back(0);
		
		for (int ii=0; ii<num_total_registers; ii++){
			bases.push_back(new_bases[ii]);
			std::cout << bases[ii] << std::endl;
		}
		
	}
	
	
	
	/**
	 get the register value at position
	 \param reggie the index of the register to get.
	 */
	int reg_val(int reggie){
		return registers[reggie];
	}
	
	/**
	 get the active register value at position
	 \param reggie the index of the ACTIVE register to get.
	 */
	int act_reg(int reggie){
		return active_registers[reggie];
	}
	
	
	/**
	 get the inactive register value at position
	 \param reggie the index of the INACTIVE register to get.
	 */
	int inact_reg(int reggie){
		return inactive_registers[reggie];
	}
	
	
	/**
	 a wrapper for incrementing.
	 
	 \return -1 if completely done. 0 if didn't roll over.  1 if rolled over, but not done, just moving to next active register.
	 
	 
	 */
	int increment(){
		
		if (DoubleOdometer::increment_registers()!=0) {
			if (DoubleOdometer::increment_active_registers()!=0)
				return -1;
			else
				return 1;
		}
		else
			return 0;
	};
	
	
	/**
	 print to cout
	 */
	void print() const{
		std::cout << "active: ";
		for (int ii=0; ii<num_active_registers; ii++)
			std::cout << active_registers[ii] << " ";
		std::cout << "\t|\t";
		
		
		std::cout << "inactive: ";
		for (int ii=0; ii<num_inactive_registers; ii++)
			std::cout << inactive_registers[ii] << " ";
		std::cout << "\t|\t";
		
		std::cout << "register values: ";
		for (int ii=0; ii<num_active_registers; ii++)
			std::cout << registers[ii] << " ";
		std::cout << "\n";
	}
};







