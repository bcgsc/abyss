#ifndef DICTIONARY_H
#define DICTIONARY_H 1

#include "ConstString.h"
#include "UnorderedMap.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

/** A bidirectional map of indices and names. */
class Dictionary
{
	public:
		typedef unsigned index_type;
		typedef unsigned index_reference;
		typedef std::string name_type;
		typedef cstring name_reference;

		typedef std::vector<const_string> Vector;
		typedef unordered_map<name_reference, index_type,
			hash<name_reference> > Map;

		Dictionary() : m_locked(false) { }

		/** Insert the specified name. */
		index_reference insert(const name_type& name)
		{
			m_vec.push_back(name);
			std::pair<Map::const_iterator, bool> inserted
				= m_map.insert(Map::value_type(
							m_vec.back(), m_map.size()));
			if (!inserted.second) {
				std::cerr << "error: duplicate ID: `"
					<< name << "'\n";
				abort();
			}
			return inserted.first->second;
		}

		/** If the specified index is within this dictionary, ensure
		 * that the name is identical, otherwise append the name to
		 * this dictionary.
		 */
		void put(index_type index, const name_type& name)
		{
			if (index < m_vec.size()) {
				assert(getName(index) == name);
			} else {
				assert(!m_locked);
				assert(index == m_vec.size());
				index_type i = insert(name);
				assert(i == index);
				(void)i;
			}
		}

		/** Return the index of the specified name. */
		index_reference getIndex(const name_type& name) const
		{
			Map::const_iterator it = m_map.find(name);
			if (it == m_map.end()) {
				std::cerr << "error: unexpected ID: `"
					<< name << "'\n";
				abort();
			}
			return it->second;
		}

		/** Return the name of the specified index. */
		name_reference getName(index_type index) const
		{
			assert(index < m_vec.size());
			return m_vec[index];
		}

		/** Lock this dictionary. No further elements may be added. */
		void lock() { m_locked = true; }

		/** Unlock this dictionary. */
		void unlock() { m_locked = false; }

		/** Return true if this dictionary is empty. */
		bool empty() const { return m_vec.empty(); }

		/** Return the number of elements in this dictionary. */
		size_t size() const { return m_vec.size(); }

		/** Return the number of elements with the specified name. */
		size_t count(const name_type& name) const
		{
			return m_map.count(name);
		}

		/** Return the last name in this dictionary. */
		name_reference back() const
		{
			assert(!m_vec.empty());
			return m_vec.back();
		}

	private:
		Map m_map;
		Vector m_vec;
		bool m_locked;
};

static inline Dictionary::name_reference get(
		const Dictionary& pmap, Dictionary::index_type index)
{
	return pmap.getName(index);
}

static inline void put(Dictionary& pmap, Dictionary::index_type index,
		const Dictionary::name_type& name)
{
	pmap.put(index, name);
}

static inline Dictionary::index_reference get(
		const Dictionary& pmap, Dictionary::name_type name)
{
	return pmap.getIndex(name);
}

#endif
