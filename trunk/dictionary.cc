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
*	FILE dictionary.cc: Implements a dictionary or associative array for pairs <string, int> .
**/

/*std includes*/
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>

/*my includes*/
#include "dictionary.h"

/*Namespaces*/
using namespace dict;

/*
* 	Static variables
* */

unsigned int Dictionary::min_capacity = 25000;

/**
 *	Dictionary function implementations
 * */

/**
 *	CONSTRUCTOR
**/
Dictionary::Dictionary()
{
	current_id = 1;
	reverse_map.resize(min_capacity);
	capacity = min_capacity;
}

/**
 *	DESTRUCTOR
**/
Dictionary::~Dictionary()
{
	table.clear();
	reverse_map.clear();
}

/**
*	FUNCTION size: Returns the size of the dictionary.  
* */
const unsigned int Dictionary::size() const
{
	return current_id - 1;
}

/**
*	FUNCTION get_term_id: Returns the id of a string or 0 in case such string is not in the dictionary.
* */

const unsigned int Dictionary::get_term_id(std::string _term) 
{
	 std::map<std::string, unsigned int, eqstr>::iterator it;
	 unsigned int id;

	 if((it = table.find(_term)) != table.end())
	 {
	 	id = it->second;
		return id;
	 }
	 else
	 {
	 	return 0;
	 }
}

/**
*	FUNCTION get_term: Returns the string associated to an id.
* */

const std::string Dictionary::get_term(unsigned int id) 
{
	return reverse_map[id-1];
}

/**
*	FUNCTION insert_term: Inserts a new string into the dictionary and returns its id.  
* */

const unsigned int Dictionary::insert_term(std::string _term)
{
	//Resizing
	if(current_id >= capacity)
	{
		reverse_map.resize(capacity + min_capacity);
		capacity += min_capacity;
	}
	
	table[_term] = current_id;

	reverse_map[current_id-1] = _term;
	current_id++;
	
	return current_id - 1;
}

