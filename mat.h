#ifndef MAT_H_
#define MAT_H_

#ifndef MAT_H_IMPLEMENTATION
#define MAT_H_IMPLEMENTATION
#endif

#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#ifndef MAT_MALLOC
#include <stdlib.h>
#define MAT_MALLOC malloc
#endif // MAT_MALLOC

#ifndef MAT_ASSERT
#include <assert.h>
#define MAT_ASSERT assert
#endif // MAT_ASSERT

typedef struct {
    size_t rows;
    size_t cols;
    size_t stride;
    double *es;
} Mat;
#define MAT_AT(m, i, j) (m).es[(i)*(m).stride + (j)]

Mat   mat_alloc(size_t rows, size_t cols);
void  mat_copy(Mat dst, Mat src);
Mat   mat_def(size_t rows, size_t cols, double *es);
void  mat_unit(Mat dst);
void  mat_fill(Mat dst, double x);
void  mat_times(Mat dst, double x);
void  mat_sum(Mat dst, Mat a);
void  mat_sub(Mat dst, Mat a);
void  mat_dot(Mat dst, Mat a, Mat b);
void  mat_print(Mat dst, const char *name, size_t padding);

#define MAT_PRINT(m) mat_print(m, #m, 0)
#define ARRAY_LEN(xs) sizeof((xs))/sizeof((xs)[0])

#endif // MAT_H_

#ifdef MAT_H_IMPLEMENTATION

Mat mat_alloc(size_t rows, size_t cols)
{
    Mat m;
    m.rows = rows;
    m.cols = cols;
    m.stride = cols;
    m.es = MAT_MALLOC(rows*cols*sizeof(*m.es));
    MAT_ASSERT(m.es != NULL);
    return m;
}

void mat_copy(Mat dst, Mat src)
{
    MAT_ASSERT(dst.rows == src.rows);
    MAT_ASSERT(dst.cols == src.cols);
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            MAT_AT(dst, i, j) = MAT_AT(src, i, j);
        }
    }
}

Mat mat_def(size_t rows, size_t cols, double *es)
{
    MAT_ASSERT(es != NULL);
    Mat m = mat_alloc(rows, cols);
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            MAT_AT(m, i, j) = es[i*rows + j];
        }
    }
    return m;
}

void mat_unit(Mat dst)
{
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            if (i == j) {
                MAT_AT(dst, i, j) = 1.f;
            }else {
                MAT_AT(dst, i, j) = 0.f;
            }
        }
    }
}

void mat_fill(Mat dst, double x)
{
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            MAT_AT(dst, i, j) = x;
        }
    }
}

void mat_times(Mat dst, double x) {
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            MAT_AT(dst, i, j) *= x;
        }
    }
}

void mat_sum(Mat dst, Mat a)
{
    MAT_ASSERT(dst.rows == a.rows);
    MAT_ASSERT(dst.cols == a.cols);
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            MAT_AT(dst, i, j) += MAT_AT(a, i, j);
        }
    }
}

void mat_sub(Mat dst, Mat a)
{
    MAT_ASSERT(dst.rows == a.rows);
    MAT_ASSERT(dst.cols == a.cols);
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            MAT_AT(dst, i, j) -= MAT_AT(a, i, j);
        }
    }
}

void mat_dot(Mat dst, Mat a, Mat b)
{
    MAT_ASSERT(a.cols == b.rows);
    MAT_ASSERT(dst.rows == a.rows);
    MAT_ASSERT(dst.cols == b.cols);
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            MAT_AT(dst, i, j) = 0;
        }
    }
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            for (size_t k = 0; k < a.cols; ++k) {
                // (i * k)  (k * j) => (i * j)
                //      ^    ^
                MAT_AT(dst, i, j) += MAT_AT(a, i, k) * MAT_AT(b, k, j);
                //                               ^              ^
            }
        }
    }
}

void mat_transpose(Mat dst, Mat a) {
    MAT_ASSERT(dst.rows == a.cols);
    MAT_ASSERT(dst.cols == a.rows);
    for (size_t i = 0; i < a.rows; i++) {
        for (size_t j = 0; j < a.cols; j++) {
            MAT_AT(dst, j, i) = MAT_AT(a, i, j);
        }
    }
}

void mat_print(Mat dst, const char *name, size_t padding)
{
    printf("%*s%s = [\n", (int) padding, "", name);
    for (size_t i = 0; i < dst.rows; ++i) {
        printf("%*s    ", (int)padding, "");
        for (size_t j = 0; j < dst.cols; ++j) {
            printf("%f ", MAT_AT(dst, i, j));
        }
        printf("\n");
    }
    printf("%*s]\n", (int) padding, "");
}


#endif // MAT_H_IMPLEMENTATION
