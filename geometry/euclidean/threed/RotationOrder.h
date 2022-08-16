#pragma once
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

 /*
  * This is not the original file distributed by the Apache Software Foundation
  * It has been modified by the Hipparchus project
  */

  //package org.hipparchus.geometry.euclidean.threed;

  //import java.util.Hash_Map;
  //import java.util.Map;

  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.geometry.Localized_Geometry_Formats;

  /**
   * This class is a utility representing a rotation order specification
   * for Cardan or Euler angles specification.
   *
   * This class cannot be instanciated by the user. He can only use one
   * of the twelve predefined supported orders as an argument to either
   * the {@link Rotation#Rotation(Rotation_Order, Rotation_Convention, double, double, double)}
   * constructor or the {@link Rotation#get_angles} method.
   *
   * sin_ce Hipparchus 1.7 this class is an enumerate class.
   *
   */
enum Rotation_Order
{
	/** Set of Cardan angles.
	 * this ordered set of rotations is around X, then around Y, then
	 * around Z
	 */
	XYZ("XYZ", Vector_3D.PLUS_I, Vector_3D.PLUS_J, Vector_3D.PLUS_K),
	/** Set of Cardan angles.
	 * this ordered set of rotations is around X, then around Z, then
	 * around Y
	 */
	 XZY("XZY", Vector_3D.PLUS_I, Vector_3D.PLUS_K, Vector_3D.PLUS_J),
	 /** Set of Cardan angles.
	  * this ordered set of rotations is around Y, then around X, then
	  * around Z
	  */
	  YXZ("YXZ", Vector_3D.PLUS_J, Vector_3D.PLUS_I, Vector_3D.PLUS_K),
	  /** Set of Cardan angles.
	   * this ordered set of rotations is around Y, then around Z, then
	   * around X
	   */
	   YZX("YZX", Vector_3D.PLUS_J, Vector_3D.PLUS_K, Vector_3D.PLUS_I),
	   /** Set of Cardan angles.
		* this ordered set of rotations is around Z, then around X, then
		* around Y
		*/
		ZXY("ZXY", Vector_3D.PLUS_K, Vector_3D.PLUS_I, Vector_3D.PLUS_J),
		/** Set of Cardan angles.
		 * this ordered set of rotations is around Z, then around Y, then
		 * around X
		 */
		 ZYX("ZYX", Vector_3D.PLUS_K, Vector_3D.PLUS_J, Vector_3D.PLUS_I),
		 /** Set of Euler angles.
		  * this ordered set of rotations is around X, then around Y, then
		  * around X
		  */
		  XYX("XYX", Vector_3D.PLUS_I, Vector_3D.PLUS_J, Vector_3D.PLUS_I),
		  /** Set of Euler angles.
		   * this ordered set of rotations is around X, then around Z, then
		   * around X
		   */
		   XZX("XZX", Vector_3D.PLUS_I, Vector_3D.PLUS_K, Vector_3D.PLUS_I),
		   /** Set of Euler angles.
			* this ordered set of rotations is around Y, then around X, then
			* around Y
			*/
			YXY("YXY", Vector_3D.PLUS_J, Vector_3D.PLUS_I, Vector_3D.PLUS_J),
			/** Set of Euler angles.
			 * this ordered set of rotations is around Y, then around Z, then
			 * around Y
			 */
			 YZY("YZY", Vector_3D.PLUS_J, Vector_3D.PLUS_K, Vector_3D.PLUS_J),
			 /** Set of Euler angles.
			  * this ordered set of rotations is around Z, then around X, then
			  * around Z
			  */
			  ZXZ("ZXZ", Vector_3D.PLUS_K, Vector_3D.PLUS_I, Vector_3D.PLUS_K),
			  /** Set of Euler angles.
			   * this ordered set of rotations is around Z, then around Y, then
			   * around Z
			   */
			   ZYZ("ZYZ", Vector_3D.PLUS_K, Vector_3D.PLUS_J, Vector_3D.PLUS_K);

/** Codes map. */
private static const Map<std::string, Rotation_Order> CODES_MAP = Hash_Map<>();
static
{
	for (const Rotation_Order type : values())
	{
		CODES_MAP.put(type.to_string(), type);
	}
}

/** Name of the rotations order. */
private const std::string name;

/** Axis of the first rotation. */
private const Vector_3D a1;

/** Axis of the second rotation. */
private const Vector_3D a2;

/** Axis of the third rotation. */
private const Vector_3D a3;

/** Private constructor.
 * This is a utility class that cannot be instantiated by the user, * so its only constructor is private.
 * @param name name of the rotation order
 * @param a1 axis of the first rotation
 * @param a2 axis of the second rotation
 * @param a3 axis of the third rotation
 */
Rotation_Order(const std::string name, const Vector_3D a1, const Vector_3D a2, const Vector_3D a3)
{
	this.name = name;
	this.a1 = a1;
	this.a2 = a2;
	this.a3 = a3;
}

/** Get a string representation of the instance.
 * @return a string representation of the instance (in fact, its name)
 */
 //override
public std::string to_string() const
{
	return name;
}

/** Get the axis of the first rotation.
 * @return axis of the first rotation
 */
public Vector_3D get_a1()
{
	return a1;
}

/** Get the axis of the second rotation.
 * @return axis of the second rotation
 */
public Vector_3D get_a2()
{
	return a2;
}

/** Get the axis of the second rotation.
 * @return axis of the second rotation
 */
public Vector_3D get_a3()
{
	return a3;
}

/**
 * Get the rotation order corresponding to a string representation.
 * @param value name
 * @return a rotation order object
 * @since 1.7
 */
public static Rotation_Order get_rotation_order(const std::string value)
{
	const Rotation_Order type = CODES_MAP.get(value);
	if (type == NULL)
	{
		// Invalid value. An exception is thrown
		throw Math_Illegal_State_Exception(Localized_Geometry_Formats.INVALID_ROTATION_ORDER_NAME, value);
	}
	return type;
}
}
