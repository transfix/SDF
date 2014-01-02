// bufferedio.cpp: implementation of the BufferedIO class.
//
//////////////////////////////////////////////////////////////////////

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "bufferedio.h"

const int BufferedIO::NDB = 1;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BufferedIO::BufferedIO(const char* fname, DiskIO::Type mode, int buf_size)
: DiskIO(mode)
{
	m_fname = STRDUP(fname);
	m_fp = NULL;
	m_ndb = buf_size;
	init();
}

BufferedIO::BufferedIO(FILE* fp, DiskIO::Type mode, int buf_size)
: DiskIO(mode)
{
	m_fname = NULL;
	m_fp = fp;
	m_ndb = buf_size;
	init();
}

void BufferedIO::init() {
	m_total = m_ndb*DBSIZE;
	m_buffer = new char[m_total];
	switch(m_mode) {
	case DiskIO::READ:
		m_bufptr = m_total;					// read pointer to buffer end
		m_bufsize = m_total;
		break;
	case DiskIO::WRITE:						// write pointer to buffer beginning
		m_bufptr = 0;						
		break;
	}
	memset(m_buffer, 0, m_total);
	m_eof = false;
}

BufferedIO::~BufferedIO()
{
	close();
	if(m_fname != NULL) free(m_fname);
	delete[] m_buffer;
}

//////////////////////////////////////////////////////////////////////
// Get Functions
//////////////////////////////////////////////////////////////////////

int BufferedIO::get(char* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	// no problem of endian
	return getraw(data, sizeof(char), n);	
}

int BufferedIO::get(unsigned char* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	// no problem of endian
	return getraw(data, sizeof(unsigned char), n);	
}

int BufferedIO::get(short* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	
#ifdef _LITTLE_ENDIAN
	short tmp;
	char* p_tmp = (char*) &tmp;
	char* p_data = (char*)data;
	n = getraw(data, sizeof(short), n);
	for(int i = 0; i < n; i++) {
		// swap byte order
		p_tmp[0] = p_data[2*i];
		p_data[2*i] = p_data[2*i+1];
		p_data[2*i+1] = p_tmp[0];
	}
#else
	n = getraw(data, sizeof(short), n);
#endif		// LITTLE_ENDIAN
	return n;
}

int BufferedIO::get(unsigned short* data, int n)
{
	return get((short*)data, n);
}

int BufferedIO::get(int* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	n = getraw(data, sizeof(int), n);
#ifdef _LITTLE_ENDIAN
	int tmp;
	char *p_tmp = (char*)&tmp;
	char *p_data = (char*)data;
	int tsize = sizeof(int);
	for(int i = 0; i < n; i++) {
		tmp = data[i];
		// swap byte order
		for(int j = tsize-1; j >= 0; j--) {
			p_data[i*tsize + tsize-j-1] = p_tmp[j];
		}
	}
#endif		// LITTLE_ENDIAN
	return n;
}

int BufferedIO::get(long* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	n = getraw(data, sizeof(long), n);
#ifdef _LITTLE_ENDIAN
	long tmp;
	char *p_tmp = (char*)&tmp;
	char *p_data = (char*)data;
	int tsize = sizeof(long);
	for(int i = 0; i < n; i++) {
		tmp = data[i];
		// swap byte order
		for(int j = tsize-1; j >= 0; j--) {
			p_data[i*tsize + tsize-j-1] = p_tmp[j];
		}
	}
#endif		// LITTLE_ENDIAN
	return n;
}

int BufferedIO::get(int64* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	n = getraw(data, sizeof(int64), n);
#ifdef _LITTLE_ENDIAN
	int64 tmp;
	char *p_tmp = (char*)&tmp;
	char *p_data = (char*)data;
	int tsize = sizeof(int64);
	for(int i = 0; i < n; i++) {
		tmp = data[i];
		// swap byte order
		for(int j = tsize-1; j >= 0; j--) {
			p_data[i*tsize + tsize-j-1] = p_tmp[j];
		}
	}
#endif		// LITTLE_ENDIAN
	return n;
}

int BufferedIO::get(float* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	n = getraw(data, sizeof(float), n);
#ifdef _LITTLE_ENDIAN
	float tmp;
	char *p_tmp = (char*)&tmp;
	char *p_data = (char*)data;
	int tsize = sizeof(float);
	for(int i = 0; i < n; i++) {
		tmp = data[i];
		// swap byte order
		for(int j = tsize-1; j >= 0; j--) {
			p_data[i*tsize + tsize-j-1] = p_tmp[j];
		}
	}
#endif		// LITTLE_ENDIAN
	return n;
}

int BufferedIO::get(double* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	//all machines have the same representation?
	return getraw(data, sizeof(double), n);	
}

//////////////////////////////////////////////////////////////////////
// Put Functions
//////////////////////////////////////////////////////////////////////

void BufferedIO::put(const char* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	putraw(data, n);
}

void BufferedIO::put(const unsigned char* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	putraw(data, n);
}

void BufferedIO::put(const short* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
#ifdef _LITTLE_ENDIAN
	short tmp;
	char* p_tmp = (char*) &tmp;
	char* p_data = (char*) data;
	int tsize = sizeof(short);
	// swap byte order
	for(int i = 0; i < n; i++) {
		for(int j = tsize-1; j >= 0; j--) {
			p_tmp[j] = p_data[i*tsize + tsize-j-1];
		}
		putraw(p_tmp, tsize);
	}
#else
	putraw(data, sizeof(short)*n);
#endif		// LITTLE_ENDIAN
}

void BufferedIO::put(const unsigned short* data, int n)
{
	put((const short*)data, n);
}

void BufferedIO::put(const int* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
#ifdef _LITTLE_ENDIAN
	int tmp;
	char* p_tmp = (char*) &tmp;
	char* p_data = (char*) data;
	int tsize = sizeof(int);
	// swap byte order
	for(int i = 0; i < n; i++) {
		for(int j = tsize-1; j >= 0; j--) {
			p_tmp[j] = p_data[i*tsize + tsize-j-1];
		}
		putraw(p_tmp, tsize);
	}
#else
	putraw(data, sizeof(int)*n);
#endif		// LITTLE_ENDIAN
}

void BufferedIO::put(const long* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
#ifdef _LITTLE_ENDIAN
	long tmp;
	char* p_tmp = (char*) &tmp;
	char* p_data = (char*) data;
	int tsize = sizeof(long);
	// swap byte order
	for(int i = 0; i < n; i++) {
		for(int j = tsize-1; j >= 0; j--) {
			p_tmp[j] = p_data[i*tsize + tsize-j-1];
		}
		putraw(p_tmp, tsize);
	}
#else
	putraw(data, sizeof(long)*n);
#endif		// LITTLE_ENDIAN
}

void BufferedIO::put(const int64* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
#ifdef _LITTLE_ENDIAN
	int64 tmp;
	char* p_tmp = (char*) &tmp;
	char* p_data = (char*) data;
	int tsize = sizeof(int64);
	// swap byte order
	for(int i = 0; i < n; i++) {
		for(int j = tsize-1; j >= 0; j--) {
			p_tmp[j] = p_data[i*tsize + tsize-j-1];
		}
		putraw(p_tmp, tsize);
	}
#else
	putraw(data, sizeof(int64)*n);
#endif		// LITTLE_ENDIAN
}

void BufferedIO::put(const float* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
#ifdef _LITTLE_ENDIAN
	float tmp;
	char* p_tmp = (char*) &tmp;
	char* p_data = (char*) data;
	int tsize = sizeof(float);
	// swap byte order
	for(int i = 0; i < n; i++) {
		for(int j = tsize-1; j >= 0; j--) {
			p_tmp[j] = p_data[i*tsize + tsize-j-1];
		}
		putraw(p_tmp, tsize);
	}
#else
	putraw(data, sizeof(float)*n);
#endif		// LITTLE_ENDIAN
}

void BufferedIO::put(const double* data, int n)
{
#ifdef _DEBUG
	assert(m_fp);
#endif
	// all machines have the same order ?!
	putraw(data, sizeof(double)*n);
}

//////////////////////////////////////////////////////////////////////
// Open/Close
//////////////////////////////////////////////////////////////////////

bool BufferedIO::open()
{
	// cannot open more than once
	if(m_fp != NULL) return true;
	switch(m_mode) {
	case DiskIO::READ:
		if((m_fp = fopen(m_fname, "r+b")) == NULL) return false;
		break;
	case DiskIO::WRITE:
		if((m_fp = fopen(m_fname, "w+b")) == NULL) return false;
		break;
	default:
		return false;
		break;
	}
	return true;
}

bool BufferedIO::reopen(const char* fname)
{
	flush();
	init();
	if(m_fp != NULL) fclose(m_fp);
	m_fp = NULL;
	if(m_fname != NULL) free(m_fname);
	m_fname = STRDUP(fname);
	return	open();
}

bool BufferedIO::reopen()
{
	if(m_fp != NULL) return true;
	switch(m_mode) {
	case DiskIO::READ:
		if((m_fp = fopen(m_fname, "r+b")) == NULL) return false;
		break;
	case DiskIO::WRITE:
		if((m_fp = fopen(m_fname, "a+b")) == NULL) return false;
		break;
	default:
		return false;
	}
	return true;
}

bool BufferedIO::close(bool fill)
{
	// flush out buffer in the write mode
	if(m_fp != NULL) {
		if(m_mode == DiskIO::WRITE) flush(fill);
		fclose(m_fp);
	}
	m_fp = NULL;
	return true;
}

void BufferedIO::flush(bool fill)
{
	switch(m_mode) {
	case DiskIO::READ:
		// !!! read header only moves to the next
		// disk block boundary
		if(m_bufptr % DBSIZE != 0) {
			m_bufptr += (DBSIZE - (m_bufptr%DBSIZE));
		}
		break;
	case DiskIO::WRITE:
		// It writes out all modified bytes up to a 
		// disk block boundary
		// only make sense for write
		if(fill) {
			if(m_bufptr % DBSIZE != 0) {
				m_bufptr += (DBSIZE - (m_bufptr%DBSIZE));
			}
		}
		if(m_bufptr != 0) fwrite(m_buffer, 1, m_bufptr, m_fp);
		// It's a good idea to zero all bytes in the buffer
		memset(m_buffer, 0, m_total);
		m_bufptr = 0;
		break;
	}
}

bool BufferedIO::seek(long offset, DiskIO::SeekMode smode)
{
	  bool sk_good; 
      if (mode() == DiskIO::READ) {
		switch(smode) {
		case DiskIO::FROM_HERE:
			// move file pointer to the real position plus the offset
		    // we optimize the seek operation if the new position is already
		    // in memory buffer
		    if((m_bufptr+offset) < m_bufsize && (m_bufptr+offset) > 0) {
			  m_bufptr += offset;
			  sk_good = true;
			} else {
			  sk_good = (fseek(m_fp, m_bufptr-m_bufsize+offset, SEEK_CUR) == 0);
			  m_bufptr = m_bufsize = m_total;
			}
			break;
		case DiskIO::FROM_START:
			sk_good = (fseek(m_fp, offset, SEEK_SET) == 0);
			// Invalidate all the data in the buffer
			m_bufptr = m_bufsize = m_total;
			break;
		case DiskIO::FROM_END:
			sk_good = (fseek(m_fp, offset, SEEK_END) == 0);
			m_bufptr = m_bufsize = m_total;
			break;
		}	    
		m_eof = false;
		return sk_good;
	  } else if (mode() == DiskIO::WRITE) {
	  flush(false);
	  switch(smode) {
	  case DiskIO::FROM_HERE:
	    return (fseek(m_fp, offset, SEEK_CUR) == 0);
	    break;
	  case DiskIO::FROM_START:
	    return (fseek(m_fp, offset, SEEK_SET) == 0);
	    break;
	  case DiskIO::FROM_END:
	    return (fseek(m_fp, offset, SEEK_END) == 0);
	    break;
	  }
	}
	return false;            // should not reach here
}

bool BufferedIO::setMode(DiskIO::Type type)
{
	if(type == mode()) return true;

	switch(type) {
	case DiskIO::READ:
		flush(false);
		m_bufptr = m_total;
		m_bufsize = m_total;
		break;
	case DiskIO::WRITE:
		fseek(m_fp, m_bufptr-m_bufsize, SEEK_CUR);
		m_bufptr = 0;
		break;
	}
	m_eof = false;
	m_mode = type;
	return true;
}

//////////////////////////////////////////////////////////////////////
// Private Helper Functions
//////////////////////////////////////////////////////////////////////

int BufferedIO::getraw(void* data, int size, int nu)
{
	if(eof()) return 0;
	// check if data already is in buffer
	if(m_bufptr + size*nu <= m_bufsize) {
		memcpy(data, m_buffer+m_bufptr, size*nu);
		m_bufptr += size*nu;
		return nu;
	}

	// return what we have if eof has been reached
	if(m_eof) {
		int n = (m_bufsize - m_bufptr) / size;
		memcpy(data, m_buffer+m_bufptr, size*n);
		m_bufptr += n*size;
		return n;
	}

	// 
	int cnt = 0;
	while(!m_eof && (m_bufptr+size * nu >= m_bufsize)){
		// add the rest of buffer to data
		int n = (m_bufsize - m_bufptr) / size;
		int re = (m_bufsize - m_bufptr) % size;   // residual bytes in buffer
		memcpy(data, m_buffer+m_bufptr, size*n);
		cnt += n;
		data = (char*)data + n*size;
		nu -= n;
		// copy the residual bytes to the front
		for(int i = 0; i < re; i++) {
			m_buffer[i] = m_buffer[m_bufsize-re+i];
		}
		m_bufptr = 0;
		// read a new buffer from disk
		if((m_bufsize = fread(m_buffer+re, 1, m_total-re, m_fp)) < m_total-re) {
			m_eof = true;
		};
		m_bufsize += re;
	}

	// recursively fill the rest of data
	if(nu) {
		return cnt + getraw(data, size, nu);
	}
	return cnt;
}

void BufferedIO::putraw(const void* data, int size)
{
	// add data to buffer if buffer is not full
	if(m_bufptr + size <= m_total) {
		memcpy(m_buffer+m_bufptr, data, size);
		m_bufptr += size;
		return;
	}

	//
	do {
		int n = m_total - m_bufptr;
		memcpy(m_buffer+m_bufptr, data, n);
		data = (char*)data + n;
		size -= n;
		m_bufptr += n;
		flush();
	} while (size >= m_total);

	if(size) {
		putraw(data, size);
	}
}



