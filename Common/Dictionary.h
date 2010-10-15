#ifndef DICTIONARY_H
#define DICTIONARY_H 1

#include "HashMap.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

/** A dictionary of strings, which gives each string a unique serial
 * number. This class serves a similar purpose to Java's
 * String::intern method.
 */
class Dictionary {
	public:
		typedef std::string key_type;
		typedef const key_type& key_reference;
		typedef std::vector<key_type> Vector;
		typedef unsigned serial_type;
		typedef hash_map<key_type, serial_type> Map;

		Dictionary() : m_locked(false) { }

		/** Insert the specified key. */
		serial_type insert(const key_type& key)
		{
			std::pair<Map::const_iterator, bool> inserted
				= m_map.insert(std::make_pair(key, m_map.size()));
			assert(inserted.second);
			m_vec.push_back(key);
			return inserted.first->second;
		}

		/** Convert the specified key to a serial number. */
		serial_type serial(const key_type& key)
		{
			std::pair<Map::const_iterator, bool> inserted
				= m_map.insert(std::make_pair(key, m_map.size()));
			if (inserted.second) {
				if (m_locked) {
					std::cerr << "error: unexpected ID: `"
						<< key << "'\n";
					exit(EXIT_FAILURE);
				}
				m_vec.push_back(key);
			}
			return inserted.first->second;
		}

		/** Convert the specified serial number to a key. */
		key_reference key(serial_type serial)
		{
			assert(serial < m_vec.size());
			return m_vec[serial];
		}

		/** Return a canonical representation for the key. */
		key_reference intern(const key_type& k)
		{
			return key(serial(k));
		}

		/** Lock this dictionary. No further keys may be added. */
		void lock() { m_locked = true; }

		/** Unlock this dictionary. */
		void unlock() { m_locked = false; }

		/** Return true if this dictionary is empty. */
		bool empty() { return m_vec.empty(); }

		/** Return the number of elements in this dictionary. */
		size_t size() { return m_vec.size(); }

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
