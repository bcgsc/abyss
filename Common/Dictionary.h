#ifndef DICTIONARY_H
#define DICTIONARY_H 1

#include "ConstString.h"
#include "HashMap.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

/** A dictionary of strings that assigns each string an index. */
class Dictionary {
	public:
		typedef std::string key_type;
		typedef cstring key_reference;
		typedef std::vector<const_string> Vector;
		typedef unsigned index_type;
		typedef hash_map<key_reference, index_type> Map;

		Dictionary() : m_locked(false) { }

		/** Insert the specified key. */
		index_type insert(const key_type& key)
		{
			m_vec.push_back(key);
			std::pair<Map::const_iterator, bool> inserted
				= m_map.insert(Map::value_type(
							m_vec.back(), m_map.size()));
			if (!inserted.second) {
				std::cerr << "error: duplicate ID: `"
					<< key << "'\n";
				exit(EXIT_FAILURE);
			}
			return inserted.first->second;
		}

		/** Return the index of the specified key. */
		index_type index(const key_type& key)
		{
			Map::const_iterator it = m_map.find(key);
			if (it != m_map.end())
				return it->second;
			if (m_locked) {
				std::cerr << "error: unexpected ID: `"
					<< key << "'\n";
				exit(EXIT_FAILURE);
			}
			return insert(key);
		}

		/** Return the name of the specified index. */
		key_reference name(index_type index)
		{
			assert(index < m_vec.size());
			return m_vec[index];
		}

		/** Lock this dictionary. No further keys may be added. */
		void lock() { m_locked = true; }

		/** Unlock this dictionary. */
		void unlock() { m_locked = false; }

		/** Return true if this dictionary is empty. */
		bool empty() { return m_vec.empty(); }

		/** Return the number of elements in this dictionary. */
		size_t size() { return m_vec.size(); }

		/** Return the number of elements with the specified key. */
		size_t count(const key_type& key) { return m_map.count(key); }

		/** Return the last key in this dictionary. */
		key_reference back()
		{
			assert(!m_vec.empty());
			return m_vec.back();
		}

	private:
		Map m_map;
		Vector m_vec;
		bool m_locked;
};

#endif
