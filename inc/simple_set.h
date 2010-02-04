#ifndef simple_set_include_guard
#define simple_set_include_guard

template<typename t>
class simple_set
{
    vector<t> elems;
    vector<int> present;
    size_t curr_size;

public:

    template<typename iter_t>
    simple_set(iter_t _begin, iter_t _end) : elems(_begin, _end)
    {
        present.resize(elems.size(), true);
        curr_size = elems.size();
    }

    simple_set()
    {
        curr_size = 0;
    }

    size_t size() const
    {
        return curr_size;
    }

    void insert(const t &val)
    {
        elems.push_back(val);
        present.push_back(1);
        curr_size++;
    }

    void revert(int index)
    {
        curr_size += !present[index];
        present[index] = 1;
    }

    void remove_index(int index)
    {
        curr_size -= !!present[index]; // !!
        present[index] = 0;
    }

	// friend class const_iterator;

    class const_iterator
    {
        const simple_set &s;
        typename vector<t>::const_iterator it;
        size_t ind;

        friend class simple_set;

        const_iterator(const simple_set &_s, size_t _ind) : s(_s), it(_s.elems.begin()), ind(_ind)
        {
            while (ind < s.present.size() && s.present[ind] == 0)
            {
                ++ind;
                ++it;
            }
        }

    public:

        bool operator!=(const const_iterator &other) const
        {
            return this->ind != other.ind;
        }

        const t& operator*() const
        {
            return *it;
        }

        const_iterator& operator++()
        {
            do
            {
                ++ind;
                ++it;
            }
            while (ind < s.present.size() && s.present[ind] == 0);

            return *this;
        }
    };

    const_iterator begin() const
    {
        return const_iterator(*this, 0);
    }

    const_iterator end() const
    {
        return const_iterator(*this, present.size());
    }
};



#endif
