/**
Copyright (c) 2011, Arlei Silva
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

@author: Arlei Silva (arleilps@gmail.com)
**/

/**
*	FILE dictionary.h: Defines a class that implements a dictionary or associative array for pairs <string, int> .
**/

#ifndef DICTIONARY_H
#define DICTIONARY_H


/*std includes*/
#include <string>
#include <vector>
#include <map>

#include <string.h>

namespace dict
{
	/**
	 *	FUNCTION cmp: Comparison function
	 * */
	struct eqstr
	{
		bool operator() (const std::string s1, const std::string s2) const
		{
			return strcmp(s1.c_str(), s2.c_str()) < 0;
		}
	};
		  
	/**
	 *	CLASS Dictionary: Maps pairs <term,id>
	 * */
	class Dictionary
	{
		public:
			/**
			 *	CONSTRUCTOR
			**/
			Dictionary();

			/**
			 *	DESTRUCTOR
			**/
			virtual ~Dictionary();
			
			/**
			*	FUNCTION getTermID: Returns the id of a string or 0 in case such string is not in the dictionary.
			* */
			const unsigned int get_term_id(std::string _term) ;
			
			/**
			*	FUNCTION getTerm: Returns the string associated to an id  .
			* */
			const std::string get_term(unsigned int id) ;
			
			/**
			*	FUNCTION insertTerm: Inserts a new string into the dictionary and returns its id.  
			* */
			const unsigned int insert_term(std::string _term);
			
			/**
			*	FUNCTION size: Returns the size of the dictionary.  
			* */
			const unsigned int size() const ;
		private:
			unsigned int current_id;
			static unsigned int min_capacity;
			unsigned int capacity;
			/*map for pairs <term, id>*/
			std::map<std::string,unsigned int,eqstr> table;
			std::vector<std::string> reverse_map; //Retrieves a term based on its id
	};
	
}

#endif

