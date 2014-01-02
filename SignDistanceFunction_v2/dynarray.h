/*
  Copyright (c): Xiaoyu Zhang (xiaoyu@csusm.edu)

  This file is part of sdf (signed distance function).

  sdf is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  sdf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with sdf; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef _DYNAMIC_ARRAY_H
#define _DYNAMIC_ARRAY_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

typedef int dyn_int;
/**
 * A simple implementation of dynamic array.
 */
template <class T>
class dynamic_array {

public:
	dynamic_array(dyn_int _size = 2);

	dynamic_array(dyn_int n, const T* _data);

	~dynamic_array();

	dyn_int length() const {
		return m_len;
	}

	dyn_int size() const {
		return m_len;
	}

	/**
	 * @return The number of elements in the array after the insertion	
	 */
	dyn_int insert(const T& x);

	dyn_int insert(const dynamic_array<T>& arry);

	dyn_int insert(dyn_int n, const T* _data);

	/*
	 *	Add n elements to the array with default value
	 */
	dyn_int increment(dyn_int n = 1);

	const T operator [] (dyn_int i) const {
		return array[i];
	}

	T& operator [] (dyn_int i);


	/**
	 * Clear all elements in the array.
	 */
	void clear();

	/**
	 * Compact the array to its minimal size.
	 * @note the minimal size is MAX(100, length())
	 */
	void compact();

	/**
	 * Expose internal data of a dynamic array.
	 * @note Use with extreme care.
	 */
	T* data() const { return array; }
private:
	T* array;
	dyn_int m_len;
	dyn_int m_size;
};

template <class T>
dynamic_array<T>::dynamic_array(dyn_int size)
{
	assert(0 <= size);
	m_size = (size < 2)? 2 : size;
	m_len = 0;

	array = (T*)malloc(sizeof(T) * m_size);
	assert(array);
}

template<class T>
dynamic_array<T>::dynamic_array(dyn_int n, const T* data)
{
	m_size = (n < 2)? 2:n;
	m_len = n;
	array = (T*)malloc(sizeof(T) * m_size);
	if(data) memcpy(array, data, sizeof(T) * m_len);
}

template <class T>
dynamic_array<T>::~dynamic_array()
{
	free(array);
}

template <class T>
dyn_int dynamic_array<T>::insert(const T& x)
{
	if(m_len == m_size) {
		m_size *= 2;
		array = (T*)realloc(array, sizeof(T)*m_size);
		assert(array);
	}
	array[m_len] = x;
	m_len++;
	return m_len;
}

template <class T>
dyn_int dynamic_array<T>::insert(dyn_int n, const T* data) 
{
	while (m_len + n > m_size) {
		m_size *= 2;
		array = (T*)realloc(array, sizeof(T)*m_size);
		assert(array);
	}
	memcpy(array+m_len, data, sizeof(T)*n);
	m_len += n;
	return m_len;
}

template <class T>
dyn_int dynamic_array<T>::increment(dyn_int n)
{
	while (m_len + n > m_size) {
		m_size *= 2;
		array = (T*)realloc(array, sizeof(T)*m_size);
		assert(array);
	}
	m_len += n;
	return m_len;
}

template<class T>
dyn_int dynamic_array<T>::insert(const dynamic_array<T>& arry) {
	return insert(arry.length(), arry.data());
}

template <class T>
T& dynamic_array<T>::operator [] (dyn_int i)
{
	if(i < m_len) {
		return array[i];
	}
	else if(i < m_size) {
		m_len = i+1;
		return array[i];
	}
	else {
		while(m_size <= i) {
			m_size = m_size << 1;
		}
		array = (T *)realloc(array, sizeof(T)*m_size);
		m_len = i+1;
	}
	return array[i];
}

template <class T>
void dynamic_array<T>::clear()
{
	m_len = 0;
}

template <class T>
void dynamic_array<T>::compact()
{
	m_size = MAX(128, length());
	array = (T *)realloc(array, sizeof(T)*m_size);
}
#endif


