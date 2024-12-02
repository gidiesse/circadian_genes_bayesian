
#include "dense_matrix.h"
#include <sstream>
#include <string>
#include <iomanip>

#include "dense_matrix.hpp"

namespace la
{
    dense_matrix::dense_matrix (size_type rows, size_type columns,
                                const_reference value)
            : m_rows (rows), m_columns (columns),
              m_data (m_rows * m_columns, value) {}

    dense_matrix::dense_matrix (std::istream & in)
    {
        read (in);
    }

    dense_matrix::size_type
    dense_matrix::sub2ind (size_type i, size_type j) const
    {
        return i * m_columns + j;
    }

    void
    dense_matrix::read (std::istream & in)
    {
        std::string line;
        std::getline (in, line);

        std::istringstream first_line (line);
        first_line >> m_rows >> m_columns;
        m_data.resize (m_rows * m_columns);

        for (size_type i = 0; i < m_rows; ++i)
        {
            std::getline (in, line);
            std::istringstream current_line (line);

            for (size_type j = 0; j < m_columns; ++j)
            {
                /* alternative syntax: current_line >> operator () (i, j);
                 * or: current_line >> m_data[sub2ind (i, j)];
                 */
                current_line >> (*this)(i, j);
            }
        }
    }

    void
    dense_matrix::swap (dense_matrix & rhs)
    {
        using std::swap;
        swap (m_rows, rhs.m_rows);
        swap (m_columns, rhs.m_columns);
        swap (m_data, rhs.m_data);
    }

    dense_matrix::reference
    dense_matrix::operator () (size_type i, size_type j)
    {
        return m_data[sub2ind (i, j)];
    }

    dense_matrix::const_reference
    dense_matrix::operator () (size_type i, size_type j) const
    {
        return m_data[sub2ind (i, j)];
    }

    dense_matrix::size_type
    dense_matrix::rows (void) const
    {
        return m_rows;
    }

    dense_matrix::size_type
    dense_matrix::columns (void) const
    {
        return m_columns;
    }

    dense_matrix
    dense_matrix::transposed (void) const
    {
        dense_matrix At (m_columns, m_rows);

        for (size_type i = 0; i < m_columns; ++i)
            for (size_type j = 0; j < m_rows; ++j)
                At(i, j) = operator () (j, i);

        return At;
    }

    dense_matrix::pointer
    dense_matrix::data (void)
    {
        return m_data.data ();
    }

    dense_matrix::const_pointer
    dense_matrix::data (void) const
    {
        return m_data.data ();
    }

    void
    dense_matrix::print (std::ostream& os) const
    {
        using size_type = dense_matrix::size_type;

        os << m_rows << " " << m_columns << "\n";

        for (size_type i = 0; i < m_rows; ++i)
        {
            for (size_type j = 0; j < m_columns; ++j)
                os << operator () (i,j) << " ";
            os << "\n";
        }
    }

    void
    dense_matrix::to_csv (std::ostream& os) const
    {
        for (size_type i = 0; i < m_rows; ++i)
        {
            for (size_type j = 0; j < m_columns; ++j)
            {
                char del = (j < m_columns - 1) ? ',' : '\n';
                os << std::setprecision(5) << operator () (i,j) << del;
            }
        }
    }

    dense_matrix
    operator * (dense_matrix const & A, dense_matrix const & B)
    {
        using size_type = dense_matrix::size_type;

        dense_matrix C (A.rows (), B.columns ());

        for (size_type i = 0; i < A.rows (); ++i)
            for (size_type j = 0; j < B.columns (); ++j)
                for (size_type k = 0; k < A.columns (); ++k)
                    C(i, j) += A(i, k) * B(k, j);

        return C;
    }

    void
    swap (dense_matrix & A, dense_matrix & B)
    {
        A.swap (B);
    }
}
