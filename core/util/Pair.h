//#pragma once
///*
// * Licensed to the Apache Software Foundation (ASF) under one or more
// * contributor license agreements.  See the NOTICE file distributed with
// * this work for additional information regarding copyright ownership.
// * The ASF licenses this file to You under the Apache License, Version 2.0
// * (the "License"); you may not use this file except in compliance with
// * the License.  You may obtain a copy of the License at
// *
// *      http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS, * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */
//
///*
// * This is not the original file distributed by the Apache Software Foundation
// * It has been modified by the Hipparchus project
// */
//
//
///**
// * Generic pair.
// * <p>
// * Although the instances of this class are immutable, it is impossible
// * to ensure that the references passed to the constructor will not be
// * modified by the caller.
// *
// * @param <K> Key type.
// * @param <V> Value type.
// */
//template<typename K, typename V>
//class Pair
//{
//private:
//    /** Key. */
//    const K key;
//    /** Value. */
//    const V value;
//
//public:
//    /**
//     * Create an entry representing a mapping from the specified key to the
//     * specified value.
//     *
//     * @param k Key (first element of the pair).
//     * @param v Value (second element of the pair).
//     */
//    Pair(K k, V v)
//    {
//        key = k;
//        value = v;
//    }
//
//    /**
//     * Create an entry representing the same mapping as the specified entry.
//     *
//     * @param entry Entry to copy.
//     */
//    Pair(Pair<? extends K, ? extends V> entry)
//    {
//        this(entry.get_key(), entry.get_value());
//    }
//
//    /**
//     * Get the key.
//     *
//     * @return the key (first element of the pair).
//     */
//    K get_key()
//    {
//        return key;
//    }
//
//    /**
//     * Get the value.
//     *
//     * @return the value (second element of the pair).
//     */
//    V get_value()
//    {
//        return value;
//    }
//
//    /**
//     * Get the first element of the pair.
//     *
//     * @return the first element of the pair.
//     */
//    K get_first()
//    {
//        return key;
//    }
//
//    /**
//     * Get the second element of the pair.
//     *
//     * @return the second element of the pair.
//     */
//    V get_second()
//    {
//        return value;
//    }
//
//    /**
//     * Compare the specified object with this entry for equality.
//     *
//     * @param o Object.
//     * @return {@code true} if the given object is also a map entry and
//     * the two entries represent the same mapping.
//     */
//    //override
//    bool equals(Object o)
//    {
//        if (this == o)
//        {
//            return true;
//        }
//		  if (dynamic_cast<const Pair*>(*o) != nullptr)
//        {
//            return false;
//        }
//        Pair<?, ?> other = (Pair<?, ?>) o;
//        return (
//            key == NULL
//                ? other.key == NULL
//                : key.equals(other.key)
//        ) && (
//            value == NULL
//                ? other.value == NULL
//                : value.equals(other.value)
//        );
//    }
//
//    /**
//     * Compute a hash code.
//     *
//     * @return the hash code value.
//     */
//    //override
//    int hash_code()
//    {
//        int result = key == NULL ? 0 : key.hash_code();
//
//        const int h = value == NULL ? 0 : value.hash_code();
//        result = 37 * result + h ^ (h >>> 16);
//
//        return result;
//    }
//
//    /** {@inherit_doc} */
//    //override
//    std::string to_string() const
//    {
//        return "[" + get_key() + ", " + get_value() + "]";
//    }
//
//    /**
//     * Convenience factory method that calls the
//     * {@link #Pair(Object, Object) constructor}.
//     *
//     * @param <K> the key type
//     * @param <V> the value type
//     * @param k First element of the pair.
//     * @param v Second element of the pair.
//     * @return a {@code Pair} containing {@code k} and {@code v}.
//     */
//    static <K, V> Pair<K, V> create(K k, V v)
//    {
//        return Pair<K, V>(k, v);
//    }
//};