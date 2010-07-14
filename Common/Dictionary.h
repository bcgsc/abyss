#ifndef DICTIONARY_H
#define DICTIONARY_H 1

#include "HashMap.h"
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
		typedef std::string Key;
		typedef unsigned Serial;
		typedef hash_map<Key, Serial> Map;

		Dictionary() : m_locked(false) { }

		/** Convert the specified key to a serial number. */
		Serial serial(const Key& key)
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
		const Key& key(Serial serial)
		{
			return m_vec.at(serial);
		}

		/** Return a canonical representation for the key. */
		const Key& intern(const Key& k)
		{
			return key(serial(k));
		}

		/** Lock this dictionary. No further keys may be added. */
		void lock() { m_locked = true; }

		/** Return true if this dictionary is empty. */
		bool empty() { return m_vec.empty(); }

		/** Return the number of elements in this dictionary. */
		size_t size() { return m_vec.size(); }

	private:
		Map m_map;
		std::vector<Key> m_vec;
		bool m_locked;
};

#endif
