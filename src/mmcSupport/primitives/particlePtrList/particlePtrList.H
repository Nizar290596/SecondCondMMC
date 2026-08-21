/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::particlePtrList
 
Description
    Container to store pointers to an object without deallocating the space 
    once the container is deconstructed.
    The OpenFOAM PtrList deallocates the stored pointers once clear() or 
    is called.


SourceFiles
    particlePtrList.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/


#ifndef mmcParticlePtrList_H
#define mmcParticlePtrList_H

#include <vector>
#include "label.H"

namespace Foam
{

template<class Type>
class particlePtrList

{
    private:
        std::vector<Type*> data_;
    
    public:
    
        //- Default constructor
        particlePtrList() = default;
        
        //- Constructor with given size
        particlePtrList(const label size)
        :   data_(size,nullptr) {}
        
        // Non-const iterator
        class iterator
        {
                //- Pointer to parent
                Type** ptr_;
            public:
            
                using iterator_category = std::forward_iterator_tag;
                using difference_type   = label;
                using value_type        = Type;
                using pointer           = value_type*;  
                using reference         = value_type&;
                friend class const_iterator;
            
                
                
                //- Constructor
                inline iterator(Type** ptr) noexcept
                : ptr_(ptr) {}
                
                inline iterator& operator++() noexcept
                {
                    ++ptr_;
                    return *this;
                }
                
                inline iterator& operator--() noexcept
                {
                    ++ptr_;
                    return *this;
                }
                
                inline pointer operator->() const
                {
                    return *ptr_;
                }
                
                inline reference operator*() const
                {
                    return **ptr_;
                }
                
                inline bool operator==(const iterator& iter) const noexcept
                {
                    return ptr_ == iter.ptr_;
                }
                inline bool operator!=(const iterator& iter) const noexcept
                {
                    return ptr_ != iter.ptr_;
                }
     
                inline bool operator<(const iterator& iter) const noexcept
                {
                    return ptr_ < iter.ptr_;
                }
                inline bool operator>(const iterator& iter) const noexcept
                {
                    return ptr_ > iter.ptr_;
                }
     
                inline bool operator<=(const iterator& iter) const noexcept
                {
                    return ptr_ <= iter.ptr_;
                }
                
                inline bool operator>=(const iterator& iter) const noexcept
                {
                    return ptr_ >= iter.ptr_;
                }
        };
        
        // Const Iterator
        class const_iterator
        {
                //- Pointer to parent
                const Type* const* ptr_;
            public: 
                
                using iterator_category = std::forward_iterator_tag;
                using difference_type   = label;
                using value_type        = const Type;
                using pointer           = const value_type*;  
                using reference         = const value_type&;
                
                //- Constructor
                inline const_iterator(Type** ptr) noexcept
                : ptr_(ptr) {}
                
                inline iterator& operator++() noexcept
                {
                    ++ptr_;
                    return *this;
                }
                
                inline iterator& operator--() noexcept
                {
                    ++ptr_;
                    return *this;
                }
                
                inline pointer operator->() const
                {
                    return *ptr_;
                }
                
                inline reference operator*() const
                {
                    return **ptr_;
                }   
        };
    
        //- Random access with automatic dererferencing
        inline Type& operator[](label i )
        {
            return *(data_[i]);
        }
        
        //- Random access with automatic dererferencing
        inline const Type& operator[](label i ) const 
        {
            return *(data_[i]);
        }
        
        //- Set the pointer at position i
        inline void set(Type* ptr, label i)
        {
            data_[i] = ptr;
        }
        
        //- Append value at end of list
        inline void append(Type* ptr)
        {
            data_.push_back(ptr);
        }
        
        //- Return size of list
        inline label size() const {return data_.size();}
        
        //- Get begin iterator
        inline iterator begin() {return iterator(&data_[0]);}
        
        //- Get end iterator
        inline iterator end() {return iterator(&data_[data_.size()]);}
        
        //- Clear the content
        void clear() {data_.clear();}
};
};  // End of namespace Foam

#endif
