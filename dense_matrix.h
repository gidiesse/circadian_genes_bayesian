//
// Created by Alessandro on 12/2/2024.
//

#ifndef BAYESIAN_STATISTICS_PROJECT_DENSE_MATRIX_H
#define BAYESIAN_STATISTICS_PROJECT_DENSE_MATRIX_H

#include <istream>
#include <vector>

namespace la // Linear Algebra
{
    class dense_matrix final
    {
        typedef std::vector<double> container_type;

    public:
        typedef container_type::value_type value_type;
        typedef container_type::size_type size_type;
        typedef container_type::pointer pointer;
        typedef container_type::const_pointer const_pointer;
        typedef container_type::reference reference;
        typedef container_type::const_reference const_reference;

    private:
        size_type m_rows, m_columns;
        container_type m_data;

        size_type
        sub2ind (size_type i, size_type j) const;

    public:
        dense_matrix (void) = default;

        dense_matrix (size_type rows, size_type columns,
                      const_reference value = 0.0);

        explicit dense_matrix (std::istream &);

        void
        read (std::istream &);

        void
        swap (dense_matrix &);

        reference
        operator () (size_type i, size_type j);
        const_reference
        operator () (size_type i, size_type j) const;

        size_type
        rows (void) const;
        size_type
        columns (void) const;

        dense_matrix
        transposed (void) const;

        pointer
        data (void);
        const_pointer
        data (void) const;

        void
        print (std::ostream& os) const;

        void
        to_csv (std::ostream& os) const;
    };

    dense_matrix
    operator * (dense_matrix const &, dense_matrix const &);

    void
    swap (dense_matrix &, dense_matrix &);

};


#endif //BAYESIAN_STATISTICS_PROJECT_DENSE_MATRIX_H
