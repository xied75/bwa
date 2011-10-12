/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <zlib.h>
#include <errno.h>
#include "utils.h"

extern time_t _prog_start;

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	// setvbuf(fp, NULL, _IOFBF, 1048576);
	return fp;
}
FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s': ", func, fn);
		perror(NULL);
		fprintf(stderr, "Abort!\n");
		abort();
	}
	// setvbuf(fp, NULL, _IOFBF, 1048576);
	return fp;
}
gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0)
		return gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
	if ((fp = gzopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	// gzbuffer(fp, 524288);
	return fp;
}
void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void err_fatal_simple_core(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}
/*
clock_t clock(void)
{
	clock_t clks = 0;
	struct timeval st;
	time_t time_now;

	// mck - use wall time ...

	gettimeofday(&st, NULL);
	time_now = st.tv_sec * 1000000L + (time_t)st.tv_usec;
	clks = (clock_t)(((double)(time_now - _prog_start) / 1000000.0) * (double)CLOCKS_PER_SEC);

	return clks;
}
*/
size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
    size_t ret = fwrite(ptr, size, nmemb, stream);
    if (ret != nmemb) 
    {
        err_fatal_simple_core("fwrite", strerror(errno));
    }
    return ret;
}

int err_printf(const char *format, ...) 
{
    va_list arg;
    int done;

    va_start(arg, format);
    done = vfprintf(stdout, format, arg);
    int saveErrno = errno;
    va_end(arg);

    if (done < 0) 
    {
        err_fatal_simple_core("vfprintf(stdout)", strerror(saveErrno));
    }
    return done;
}

int err_fprintf(FILE *stream, const char *format, ...) 
{
    va_list arg;
    int done;

    va_start(arg, format);
    done = vfprintf(stream, format, arg);
    int saveErrno = errno;
    va_end(arg);

    if (done < 0) 
    {
        err_fatal_simple_core("vfprintf", strerror(saveErrno));
    }
    return done;
}

int err_fflush(FILE *stream) 
{
    int ret = fflush(stream);
    if (ret != 0) 
    {
        err_fatal_simple_core("fflush", strerror(errno));
    }
    return ret;
}

int err_fclose(FILE *stream) 
{
    int ret = fclose(stream);
    if (ret != 0) 
    {
        err_fatal_simple_core("fclose", strerror(errno));
    }
    return ret;
}
