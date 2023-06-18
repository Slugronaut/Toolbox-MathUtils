/*
 * Copyright 2012-2015
 * James Clark
 */
using UnityEngine;
using System;
using System.Diagnostics;
using UnityEngine.Assertions;
using System.Collections.Generic;
using System.Linq;
using System.Collections;


namespace Toolbox
{
    /// <summary>
    /// Timing utility class.
    /// </summary>
    public static class Timing
    {
        public static TimeSpan Time(Action action)
        {
            Stopwatch stopwatch = Stopwatch.StartNew();
            action();
            stopwatch.Stop();
            return stopwatch.Elapsed;
        }

        public static long TimeInMilli(Action action)
        {
            Stopwatch stopwatch = Stopwatch.StartNew();
            action();
            stopwatch.Stop();
            return stopwatch.ElapsedMilliseconds;
        }
    }
}

namespace Toolbox.Math
{
    /// <summary>
    /// Common math tools for Lerping, Vector manipulation, collision detection, and physics integration.
    /// </summary>
	public class MathUtils
    {
        #region ScreenSpace
        public static Vector3 GetMouseWorldPos()
        {
            Vector3 mousePos = Input.mousePosition;
            mousePos.z = -Camera.main.transform.position.z;
            return Camera.main.ScreenToWorldPoint(mousePos);
        }
        #endregion


        #region Lerping
        public enum EaseType
        {
            EaseInOut,
            EaseOut,
            EaseIn,
            Linear
        }

        /// <summary>
        /// Converts a linear volume between 0 and 1 to a value between -80 and 0 decibles.
        /// </summary>
        /// <param name="linear"></param>
        /// <returns></returns>
        public static float LinearToDecible(float linear)
        {
            return 20 * Mathf.Log10(linear);
        }

        /// <summary>
        /// Converts a decible volume between -80 and 0 to a linear volume between 0 and 1.
        /// </summary>
        /// <param name="decible"></param>
        /// <returns></returns>
        public static float DecibleToLinear(float decible)
        {
            return Mathf.Pow(10.0f, decible / 20);
        }

        /// <summary>
        /// Convienience function that allows us to lerp based on an enumerated smoothing type.
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <param name="value"></param>
        /// <param name="type"></param>
        /// <returns></returns>
        public static float EaseFromTo(float start, float end, float value, EaseType type = EaseType.EaseInOut)
        {
            value = Mathf.Clamp01(value);

            switch (type)
            {
                case EaseType.EaseInOut:
                    return Mathf.Lerp(start, end, value * value * (3.0f - 2.0f * value));

                case EaseType.EaseOut:
                    return Mathf.Lerp(start, end, Mathf.Sin(value * Mathf.PI * 0.5f));

                case EaseType.EaseIn:
                    return Mathf.Lerp(start, end, 1.0f - Mathf.Cos(value * Mathf.PI * 0.5f));

                default:
                    return Mathf.Lerp(start, end, value);
            }
        }

        /// <summary>
        /// Manual smoothing function for a single dimension.
        /// </summary>
        /// <param name="pastPosition"></param>
        /// <param name="pastTargetPosition"></param>
        /// <param name="targetPosition"></param>
        /// <param name="speed"></param>
        /// <param name="deltaTime"></param>
        /// <returns></returns>
        public static float SmoothApproach(float pastPosition, float pastTargetPosition, float targetPosition, float speed, float deltaTime)
        {
            float t = deltaTime * speed;
            float v = (targetPosition - pastTargetPosition) / Mathf.Max(t, 0.00001f);
            float f = pastPosition - pastTargetPosition + v;
            float result = targetPosition - v + f * Mathf.Exp(-t);
            return result;
        }

        /// <summary>
        /// Scales the input and output of a curve to fit within a min and max range.
        /// </summary>
        /// <param name="curve"></param>
        /// <param name="input"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        static public float EvaluateCurveRange(float input, float min, float max, AnimationCurve curve)
        {
            float diff = max - min;
            return (curve.Evaluate(((input - min) / diff)) * diff) + min;
        }

        /// <summary>
        /// normalizes the input within the range of min and max and evaluates a curve with it.
        /// </summary>
        /// <param name="curve"></param>
        /// <param name="input"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        static public float EvaluateCurve(float input, float min, float max, AnimationCurve curve)
        {
            var diff = max - min;
            return curve.Evaluate((input - min) / diff);
        }

        /// <summary>
        /// Evaluates a value between min and max using a sampled percentage of a curve.
        /// </summary>
        /// <param name="percent"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <param name="curve"></param>
        /// <returns></returns>
        static public float EvaluateCurveScale(float percent, float min, float max, AnimationCurve curve)
        {
            var scale = curve.Evaluate(percent);
            float diff = max - min;
            return min + (diff * scale);
        }

        /// <summary>
        /// Samples a random point on an AnimationCurve.
        /// </summary>
        /// <param name="curve"></param>
        /// <returns></returns>
        public static float SampleRandomPointOnCurve(AnimationCurve curve)
        {
            var keys = curve.keys;
            var finalKey = keys[keys.Length - 1];
            return curve.Evaluate(UnityEngine.Random.Range(0, finalKey.time));
        }

        /// <summary>
        /// Samples a random point on an AnimationCurve.
        /// </summary>
        /// <param name="curve"></param>
        /// <returns></returns>
        public static float SampleRandomPointOnCurve(AnimationCurve curve, float min)
        {
            var keys = curve.keys;
            var finalKey = keys[keys.Length - 1];
            return curve.Evaluate(UnityEngine.Random.Range(Mathf.Max(min, keys[0].time), finalKey.time));
        }

        /// <summary>
        /// Samples a random point on an AnimationCurve.
        /// </summary>
        /// <param name="curve"></param>
        /// <returns></returns>
        public static float SampleRandomPointOnCurve(AnimationCurve curve, float min, float max)
        {
            var keys = curve.keys;
            var finalKey = keys[keys.Length - 1];
            return curve.Evaluate(UnityEngine.Random.Range(Mathf.Max(min, keys[0].time), Mathf.Min(max, finalKey.time)));
        }

        /// <summary>
        /// Converts an input within a min and max range to the corresponding
        /// value within the range of the output min and max range
        /// </summary>
        /// <param name="input"></param>
        /// <param name="sourceMin"></param>
        /// <param name="sourceMax"></param>
        /// <returns></returns>
        public static float ConvertRange(float input, float inMin, float inMax, float outMin, float outMax)
        {
            float outRange = outMax - outMin;
            float inNorm = input / (inMax - inMin);
            return (inNorm * outRange) + outMin;
        }

        /// <summary>
        /// Converts an input within a min and max range to the corresponding
        /// value within the range of the output max and min range. Note that the
        /// output max and min ranges are inverted so that the input min corresponds
        /// to the output max and visa versa.
        /// </summary>
        /// <param name="input"></param>
        /// <param name="sourceMin"></param>
        /// <param name="sourceMax"></param>
        /// <returns></returns>
        public static float ConvertInvertedRange(float input, float inMin, float inMax, float outMin, float outMax)
        {
            float omx = Mathf.Max(outMax, outMin);
            float omn = Mathf.Min(outMax, outMin);
            float imx = Mathf.Max(inMax, inMin);
            float imn = Mathf.Min(inMax, inMin);

            float outRange = omx - omn;
            float inNorm = input / (imx - imn);
            return (inNorm * ((outMax > outMin) ? outRange : 1.0f/outRange) ) + omn;
        }

        /// <summary>
        /// Returns a non-clamped normalized value based on an input, a min, and a max.
        /// </summary>
        /// <param name="input"></param>
        /// <param name="outMax"></param>
        /// <param name="inMax"></param>
        /// <returns></returns>
        public static float Normalize(float input, float inMin, float inMax)
        {
            Assert.IsTrue(inMax >= inMin);
            return input / (inMax - inMin);
        }

        /// <summary>
        /// Converts an input within a min and max range to the corresponding
        /// value within the range of the output min and max range as a normalized value.
        /// </summary>
        /// <param name="input"></param>
        /// <param name="sourceMin"></param>
        /// <param name="sourceMax"></param>
        /// <returns></returns>
        public static float ConvertRangeNormalized(float input, float inMin, float inMax, float outMin, float outMax)
        {
            float r = inMax - inMin;
            float norm = input / r;
            return (norm * (outMax - outMin)) + outMin;
        }

        /// <summary>
        /// Returns an unclamped value of 0 to 1 that represents the input's percentage along the range of min to max.
        /// </summary>
        /// <param name="input"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        public static float PercentageOfRange(float input, float min, float max)
        {
            Assert.IsTrue(max >= min);
            return (input - min) / (max - min);
        }

        /**
	     * Simple tool for linear interpolation.
	     *
	     * @param last
	     * @param now
	     * @param alpha
	     * @return
	     */
        private float SimpleLerp(float last, float now, float alpha)
        {
            //alpha * black + (1 - alpha) * red
            return (alpha * now) + ((1 - alpha) * last);
        }

        /**
	     * Basic linear interpolation based on elapsed time.
	     * 
	     * @param start Starting value of the lerp
	     * @param target Ending value of the lerp
	     * @param duration Length of time for the lerp (the max lerp time)
	     * @param timeSinceStart Current time since lerp started
	     * @return 
	     */
        public static float TimeLerp(float start, float target, float duration, float timeSinceStart)
        {
            float value = start;
            if (timeSinceStart > 0 && timeSinceStart < duration)
            {
                float range = target - start;
                float percent = timeSinceStart / duration;
                value = start + (range * percent);
            }
            else if (timeSinceStart >= duration)
            {
                value = target;
            }
            return value;
        }

        /**
	     * Basic spring based on elapsed time.
	     * 
	     * @param start
	     * @param target
	     * @param duration
	     * @param timeSinceStart
	     * @return 
	     */
        public static float TimeEase(float start, float target, float duration, float timeSinceStart)
        {
            float value = start;
            if (timeSinceStart > 0.0f && timeSinceStart < duration)
            {
                float range = target - start;
                float percent = timeSinceStart / (duration * 0.5f);
                if (percent < 1.0f)
                {
                    value = start + ((range * 0.5f) * percent * percent * percent);
                }
                else
                {
                    float shiftedPercent = percent - 2.0f;
                    value = start + ((range * 0.5f)
                            * ((shiftedPercent * shiftedPercent * shiftedPercent) + 2.0f));
                }
            }
            else if (timeSinceStart >= duration)
            {
                value = target;
            }
            return value;
        }

        public static Vector3 SmoothApproach(Vector3 pastPosition, Vector3 pastTargetPosition, Vector3 targetPosition, float speed)
        {
            if (!(Time.deltaTime > 0))
                return pastPosition;
            if (!(speed > 0))
                speed = float.Epsilon;

            float t = Time.deltaTime * speed;
            Vector3 v = (targetPosition - pastTargetPosition) / t;
            Vector3 f = pastPosition - pastTargetPosition + v;
            return targetPosition - v + f * Mathf.Exp(-t);
        }

        public static Vector2 SmoothApproach(Vector2 pastPosition, Vector2 pastTargetPosition, Vector2 targetPosition, float speed)
        {
            if (!(Time.deltaTime > 0))
                return pastPosition;
            if (!(speed > 0))
                speed = float.Epsilon;

            float t = Time.deltaTime * speed;
            Vector2 v = (targetPosition - pastTargetPosition) / t;
            Vector2 f = pastPosition - pastTargetPosition + v;
            return targetPosition - v + f * Mathf.Exp(-t);
        }

        public static float SmoothApproach(float pastPosition, float pastTargetPosition, float targetPosition, float speed)
        {
            if (!(Time.deltaTime > 0))
                return pastPosition;
            if (!(speed > 0))
                speed = 0.0001f;
            
            float t = Time.deltaTime * speed;
            float v = (targetPosition - pastTargetPosition) / t;
            float f = pastPosition - pastTargetPosition + v;
            return targetPosition - v + f * Mathf.Exp(-t);
        }
        #endregion


        #region Collision
        /// <summary>
        /// Is a given world-space position within a camera's viewport.
        /// </summary>
        /// <param name="cam">The camer whos viewport is being checked.</param>
        /// <param name="agentPos">The word-space position that is being queried.</param>
        /// <param name="safeX">An additional 'safe zone' buffer for the viewport's width.</param>
        /// <param name="safeY">An additional 'safe zone' buffer for the viewport's height.</param>
        /// <returns></returns>
        public static bool IsInViewport(Camera cam, Vector3 agentPos, float safeX = 0, float safeY = 0)
        {
            Vector3 screenPoint = cam.WorldToViewportPoint(agentPos);
            return screenPoint.x > safeX && screenPoint.x < 1 - safeX &&
                                screenPoint.y > safeY && screenPoint.y < 1 - safeY;
        }

        /// <summary>
        /// Is the given position in world-space occluded by geometry in respects to a given camera.
        /// </summary>
        /// <param name="cam"></param>
        /// <param name="agentPos"></param>
        /// <param name="radius"></param>
        /// <param name="mask"></param>
        /// <param name="triggerInteraction"></param>
        /// <param name="debug"></param>
        /// <returns></returns>
        public static bool IsOcculded(Camera cam, Vector3 agentPos, float radius, LayerMask mask, QueryTriggerInteraction triggerInteraction, bool debug = false)
        {
            Vector3 origin = cam.transform.position;
            Vector3 dir = agentPos - origin;
            float dist = dir.magnitude;
            dir.Normalize();

            bool result = Physics.SphereCast(origin, radius, dir, out RaycastHit info, dist, mask, triggerInteraction);

            #if UNITY_EDITOR
            if (debug)
            {
                if (result)
                    UnityEngine.Debug.DrawLine(origin, info.point, Color.red);
                else UnityEngine.Debug.DrawLine(origin, agentPos, Color.white);
                return !result;
            }
            else return !result;
            #else
            return !result;
            #endif
        }

        /// <summary>
        /// Returns a bounding box that contains everything in the object.
        /// </summary>
        /// <returns></returns>
        public static Bounds GetObjectBounds(Vector3 center, GameObject gameObject, bool includeRenders, bool includeColliders, bool includeColliders2D)
        {
            Bounds b = new Bounds(center, Vector3.zero);
            if (includeRenders)
            {
                var rends = gameObject.GetComponentsInChildren<Renderer>(true);
                foreach (Renderer r in rends)
                    b.Encapsulate(r.bounds);
            }

            if (includeColliders)
            {
                var cols = gameObject.GetComponentsInChildren<Collider>(true);
                foreach (Collider c in cols)
                    b.Encapsulate(c.bounds);
            }

            if (includeColliders2D)
            {
                var cols = gameObject.GetComponentsInChildren<Collider2D>(true);
                foreach (Collider2D c in cols)
                    b.Encapsulate(c.bounds);
            }
            return b;
        }

        /// <summary>
        /// Returns a bounding box that contains everything in the object.
        /// Only includes components that are contained within the specified LayerMask.
        /// </summary>
        /// <returns></returns>
        public static Bounds GetObjectBounds(Vector3 center, GameObject gameObject, bool includeRenders, bool includeColliders, bool includeColliders2D, LayerMask colliderLayers)
        {
            Bounds b = new Bounds(center, Vector3.zero);
            if (includeRenders)
            {
                var rends = gameObject.GetComponentsInChildren<Renderer>(false);
                foreach (Renderer r in rends)
                {
                    if(colliderLayers.ContainsLayerIndex(r.gameObject.layer))
                        b.Encapsulate(r.bounds);
                }
            }

            if (includeColliders)
            {
                var cols = gameObject.GetComponentsInChildren<Collider>(false);
                foreach (Collider c in cols)
                {
                    if(colliderLayers.ContainsLayerIndex(c.gameObject.layer))
                        b.Encapsulate(c.bounds);
                }
            }

            if (includeColliders2D)
            {
                var cols = gameObject.GetComponentsInChildren<Collider2D>(false);
                foreach (Collider2D c in cols)
                {
                    if(colliderLayers.ContainsLayerIndex(c.gameObject.layer))
                        b.Encapsulate(c.bounds);
                }
            }
            return b;
        }

        /// <summary>
        /// Returns a bounding box that contains everything in the scene.
        /// </summary>
        /// <returns></returns>
        public static Bounds GetAllBounds(bool includeRenders, bool includeColliders, bool includeColliders2D)
        {
            Bounds b = new Bounds(Vector3.zero, Vector3.zero);
            if (includeRenders)
            {
                foreach (Renderer r in GameObject.FindObjectsOfType(typeof(Renderer)))
                    b.Encapsulate(r.bounds);
            }

            if (includeColliders)
            {
                foreach (Collider c in GameObject.FindObjectsOfType(typeof(Collider)))
                    b.Encapsulate(c.bounds);
            }

            if (includeColliders2D)
            {
                foreach (Collider2D c in GameObject.FindObjectsOfType(typeof(Collider2D)))
                    b.Encapsulate(c.bounds);
            }
            return b;
        }

        /// <summary>
        /// Returns the closest collider vertex to a collision
        /// that occured on a polygon collider.
        /// </summary>
        /// <param name="hit"></param>
        /// <returns></returns>
        public static Vector3 ClosestColliderVertex(RaycastHit2D hit)
        {
            PolygonCollider2D pc = hit.collider as PolygonCollider2D;

            float minDistanceSqr = Mathf.Infinity;
            Vector3 nearestVertex = Vector3.zero;

            // Scan all collider points to find nearest
            foreach (Vector3 colliderPoint in pc.points)
            {
                // Convert to world point
                Vector3 colliderPointWorld = hit.transform.TransformPoint(colliderPoint);

                Vector3 diff = hit.point - (Vector2)colliderPointWorld;
                float distSqr = diff.sqrMagnitude;

                if (distSqr < minDistanceSqr)
                {
                    minDistanceSqr = distSqr;
                    nearestVertex = colliderPointWorld;
                }
            }

            return nearestVertex;
        }

        /// <summary>
        /// Tests to see if a line segment intersects a circle.
        /// </summary>
        /// <param name="lineStart"></param>
        /// <param name="lineEnd"></param>
        /// <param name="circleCenter"></param>
        /// <param name="circleRadius"></param>
        /// <returns></returns>
        public static bool LineIntersectsCircle(Vector2 lineStart, Vector2 lineEnd, Vector2 circleCenter, float circleRadius)
        {
            Vector2 dir = lineEnd - lineStart;
            Vector2 fuc = lineStart - circleCenter;

            float a = Vector2.Dot(dir, dir);
            float b = 2 * Vector2.Dot(fuc, dir);
            float c = Vector2.Dot(fuc, fuc) - (circleRadius * circleRadius);

            float disc = b * b - 4 * a * c;
            if (disc < 0.0f) return false;

            disc = Mathf.Sqrt(disc);
            //float t1 = (-b - disc) / (2 * a);
            //float t2 = (-b + disc) / (2 * a);
            float delt = 1 / (2 * a);
            float t1 = (-b - disc) * delt;
            float t2 = (-b + disc) * delt;

            if (t1 >= 0.0f && t1 <= 1.0f) return true;
            if (t2 >= 0.0f && t2 <= 1.0f) return true;

            return false;
        }

        /// <summary>
        /// Takes a collision that occurred within a polygon and finds the closest point along
        /// the outside edge of that collider. The collection of points is expected to be in local space.
        /// </summary>
        /// <param name="hit">The collision event that happened within a polygon collider.</param>
        /// <param name="offset">How far outside from the cloest edge we should project. Useful if you need to find a location slightly outside of the polygon edge.</param>
        /// <returns></returns>
        public static Vector3 ClosestColliderEdgePoint(Vector3 point, Vector2[] points, Transform localOrigin, float offset)
        {
            //PolygonCollider2D pc = hit.collider as PolygonCollider2D;
            Vector3 nearestPoint = Vector3.zero;
            float nearestDist = float.MaxValue;
            float dist;
            Vector3 p0, p1, p2;

            if (points == null) return point;

            for (int i = 0; i < points.Length - 1; i++)
            {
                p1 = localOrigin.TransformPoint(points[i]);
                p2 = localOrigin.TransformPoint(points[i + 1]);
                p0 = ClosestPointToLine(p1, p2, point);
                dist = (point - p0).sqrMagnitude; //Vector3.Distance(p0, point);
                dist *= dist;
                if (dist <= nearestDist)
                {
                    nearestDist = dist;
                    nearestPoint = p0;
                }
            }
            //edge case, test last line using last and first vertex
            p1 = localOrigin.TransformPoint(points[points.Length - 1]);
            p2 = localOrigin.TransformPoint(points[0]);
            p0 = ClosestPointToLine(p1, p2, point);
            dist = (point - p0).sqrMagnitude; //Vector3.Distance(p0, point);
            dist *= dist;
            if (dist <= nearestDist)
            {
                nearestDist = dist;
                nearestPoint = p0;
            }

            if (Mathf.Abs(offset) > 0.0f)
            {
                //Vector3 surface = nearestPoint;
                Vector2 proj = (Vector2)nearestPoint - (Vector2)point;
                nearestPoint = ((Vector2)nearestPoint + ((proj.normalized) * offset));
                //Debug.DrawLine(hit.point, surface, Color.green, 0.25f);
                //Debug.DrawLine(surface, nearestPoint, Color.magenta, 0.25f);
            }
            return nearestPoint;
        }

        /// <summary>
        /// Checks for line-of-sight between two points.
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <param name="mask"></param>
        /// <returns></returns>
        public static bool HasLos(Vector3 start, Vector3 end, LayerMask mask)
        {
            return !Physics.Linecast(start, end, mask, QueryTriggerInteraction.Ignore);
        }

        /// <summary>
        /// Checks for line-of-sight between two points.
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <param name="mask"></param>
        /// <returns></returns>
        public static bool HasLos(Vector3 start, Vector3 end, float radius, float maxDist, LayerMask mask)
        {
            var dir = (end - start);
            maxDist = Mathf.Min(maxDist, dir.magnitude);

            if (radius == 0)
                return !Physics.Raycast(new Ray(start, dir.normalized), maxDist, mask, QueryTriggerInteraction.Ignore);
            else return !Physics.SphereCast(new Ray(start, dir.normalized), radius, maxDist, mask, QueryTriggerInteraction.Ignore);
        }
        #endregion


        #region Vector and Angles
        /// <summary>
        /// Given a source, target, and trajectory angle, this will calculate the needed velocity to reach the desired location.
        /// It takes global gravity scale into account automatically.
        /// </summary>
        /// <param name="source"></param>
        /// <param name="target"></param>
        /// <param name="angle"></param>
        /// <returns></returns>
        public static Vector3 BallisticVelocityVector(Vector3 source, Vector3 target, float angle)
        {
            Vector3 direction = target - source;            // get target direction
            float h = direction.y;                                            // get height difference
            direction.y = 0;                                                // remove height
            float distance = direction.magnitude;                            // get horizontal distance
            float a = angle * Mathf.Deg2Rad;                                // Convert angle to radians
            direction.y = distance * Mathf.Tan(a);                            // Set direction to elevation angle
            distance += h / Mathf.Tan(a);                                        // Correction for small height differences

            // calculate velocity
            float velocity = Mathf.Sqrt(distance * Physics.gravity.magnitude / Mathf.Sin(2 * a));
            //this is a guard agaibst bullshit
            if (Mathf.Abs(velocity) < 0.01f)
                velocity = 0.01f;
            return velocity * direction.normalized;

        }

        /// <summary>
        /// Checks to see if a point is within a view angle of another point on the X,Z plane.
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        public static bool InAbsoluteAngleXZ(Vector3 agent, Vector3 target, float viewAngle)
        {
            var targetPos = new Vector2(target.x, target.z);
            var agentPos = new Vector2(agent.x, agent.z);

            float left = Vector2.Angle(targetPos - agentPos, Vector2.left);
            float right = Vector2.Angle(targetPos - agentPos, Vector2.right);
            return (left < viewAngle && left > 0) ||
                   (right < viewAngle && right > 0);
        }

        /// <summary>
        /// Returns the rotation angle on the Y-axis needed to rotate agent towards target.
        /// </summary>
        /// <param name="agent"></param>
        /// <param name="target"></param>
        /// <returns></returns>
        public static float AngleOnY(Vector3 agentPos, Vector3 agentForward, Vector3 targetPos)
        {
            var target2D = new Vector3(targetPos.x, targetPos.z);
            var agent2D = new Vector3(agentPos.x, agentPos.z);
            var targetDir = target2D - agent2D;
            return Vector2.Angle(agentForward, targetDir);
        }

        /// <summary>
        /// Gets quaternion that represents a rotation from agentPos to targetPos with an optional offset angle.
        /// </summary>
        /// <param name="targetPos"></param>
        /// <param name="agentPos"></param>
        /// <param name="angleOffset"></param>
        /// <returns></returns>
        public static Quaternion LookRotation(Vector3 agentPos, Vector3 targetPos, Vector3 angleOffset)
        {
            var forward = targetPos - agentPos;
            var rot = Quaternion.LookRotation(forward);

            if (angleOffset != Vector3.zero)
            {
                var eangles = rot.eulerAngles;
                eangles.x += angleOffset.x;
                eangles.y += angleOffset.y;
                eangles.z += angleOffset.z;
                rot.eulerAngles = eangles;
            }

            return rot;
        }

        /// <summary>
        /// Gets quaternion that represents a rotation from agentPos to targetPos with an optional offset angle.
        /// </summary>
        /// <param name="targetPos"></param>
        /// <param name="agentPos"></param>
        /// <param name="angleOffset"></param>
        /// <returns></returns>
        public static Quaternion LookRotation(Vector3 agentPos, Vector3 targetPos)
        {
            var forward = targetPos - agentPos;
            var rot = Quaternion.LookRotation(forward);
            return rot;
        }

        /// <summary>
        /// Adds a randomness factor to a directional vector based on spread.
        /// The vector is assumed to be in aim-space and will only randomize the x and y elements.
        /// Spread can be between 0 and 1. While not technically accurate, this
        /// calculation is very fast. For better accuracy, use RandomInCone.
        /// </summary>
        /// <param name="forward"></param>
        /// <param name="maxAngle"></param>
        /// <returns></returns>
        public static Vector3 RandomForwardSpreadFast(Vector3 forward, float maxAngle)
        {
            //return forward + UnityEngine.Random.insideUnitSphere * Mathf.Clamp01(spread);
            var lookRot = Quaternion.LookRotation(forward);
            var randRot = UnityEngine.Random.rotation;
            lookRot = Quaternion.RotateTowards(lookRot, randRot, UnityEngine.Random.Range(0, maxAngle * 0.5f));
            return lookRot * forward;
        }

        /// <summary>
        /// Gives a random vector within a cone that has a given radius.
        /// The vector is assumed to be in aim-space and will only randomize the x and y elements.
        /// Computationally expensive, for better speed use RandomForwardSpreadFast.
        /// </summary>
        /// <param name="radius"></param>
        /// <returns></returns>
        public static Vector3 RandomInCone(float radius)
        {
            float z = UnityEngine.Random.Range(Mathf.Cos(radius), 1);
            float t = UnityEngine.Random.Range(0, Mathf.PI * 2);
            return new Vector3(Mathf.Sqrt(1 - z * z) * Mathf.Cos(t), Mathf.Sqrt(1 - z * z) * Mathf.Sin(t), z);
        }

        /// <summary>
        /// Offsets a vector in forward-space: x-axis is always left/right,
        /// y-axis is always up/down, and z-axis is always forward/back
        /// relative to 'dir'.
        /// </summary>
        /// <param name="pos"></param>
        /// <param name="dir"></param>
        /// <param name="offset"></param>
        /// <returns></returns>
        public static Vector3 ForwardSpaceOffset(Vector3 pos, Vector3 dir, Vector3 offset)
        {
            return pos + Vector3.Cross(dir, new Vector3(0, -offset.x, 0)) +
                    new Vector3(0, offset.y, 0) +
                    (dir * offset.z);
        }

        /// <summary>
        /// Helper method for determining if a point is to the left or right side of a line.
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static bool IsLeft(Vector2 start, Vector2 end, Vector2 p)
        {
            return ((end.x - start.x) * (p.y - start.y) - (end.y - start.y) * (p.x - start.x)) > 0;
        }

        public enum Side
        {
            Center,
            Left,
            Right,
        }

        /// <summary>
        /// Faster way of calculating angle between two vectors.
        /// As a bonus it Returns a value between -180 and 180 which
        /// allows one to determine which side the 'to' angle was on
        /// without requiring an additional cross product!
        /// </summary>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <returns></returns>
        public static float FastAngle(Vector2 from, Vector2 to)
        {
            float a = Mathf.Atan2(from.y, from.x) * Mathf.Rad2Deg;
            float b = Mathf.Atan2(to.y, to.x) * Mathf.Rad2Deg;
            return Mathf.DeltaAngle(a, b);
        }

        /// <summary>
        /// Helper that can be used to determine if a point is on the right or left side of a Transform.
        /// Note that Left and Right don't necessarily mean left and right side of a Transform's forward vector.
        /// It's just an easy way to determine if we are on one side of a vector-angle comparison or the other.
        /// </summary>
        /// <param name="anchor"></param>
        /// <param name="target"></param>
        /// <returns></returns>
        public static Side GetSide(Transform anchor, Vector3 target)
        {
            var relPoint = anchor.InverseTransformPoint(target);
            if (relPoint.x < 0.0) return Side.Left;
            else if (relPoint.x > 0) return Side.Right;
            return Side.Center;
        }

        /// <summary>
        /// Rotates a 2D vector by angle.
        /// </summary>
        /// <param name="v"></param>
        /// <param name="angle"></param>
        /// <returns></returns>
        public static Vector2 Rotate(Vector2 v, float angle)
        {
            float sin = Mathf.Sin(angle);
            float cos = Mathf.Cos(angle);
            float tx = v.x;
            float ty = v.y;
            return new Vector2((cos * tx) - (sin * ty),
                               (cos * ty) + (sin * tx));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="inDirection"></param>
        /// <param name="inNormal"></param>
        /// <returns></returns>
        public static Vector2 Reflect(Vector2 inDirection, Vector2 inNormal)
        {
            return 2.0f * Vector2.Dot(inDirection, inNormal) * inNormal - inDirection;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <param name="testPoint"></param>
        /// <returns></returns>
        public static Vector3 ClosestPointToLine(Vector3 start, Vector3 end, Vector3 testPoint)
        {
            Vector3 pointTowardStart = testPoint - start;
            Vector3 startTowardEnd = (end - start).normalized;

            float lengthOfLine = Vector3.Distance(start, end);
            float dotProduct = Vector3.Dot(startTowardEnd, pointTowardStart);

            if (dotProduct <= 0)
                return start;

            if (dotProduct >= lengthOfLine)
                return end;

            Vector3 thirdVector = startTowardEnd * dotProduct;
            Vector3 closestPointOnLine = start + thirdVector;
            return closestPointOnLine;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="circleCenter"></param>
        /// <param name="arcStart"></param>
        /// <param name="arcEnd"></param>
        /// <param name="testPoint"></param>
        /// <returns></returns>
        public static Vector3 ClosestPointToArc(Vector3 circleCenter, Vector3 arcStart, Vector3 arcEnd, Vector3 testPoint)
        {
            //Vector3 endTowardBegin = _arcEnd - _arcBegin;
            float radius = (arcStart - circleCenter).magnitude;
            Vector3 surfacePoint = circleCenter - ((circleCenter - testPoint).normalized * radius);

            //if (Vector3.Dot(Perpendicular(endTowardBegin), surfacePoint - _arcBegin) > 0)
            if (Vector3.Dot(Vector3.Cross(arcEnd, arcStart).normalized, surfacePoint - arcStart) > 0)
            {
                // The closest point on the circle is contained within the arc.
                return surfacePoint;
            }
            else
            {
                // The closest point on the circle is outside the arc. One of the endpoints is closest.
                float distanceToBegin = Vector3.Distance(testPoint, arcStart);
                float distanceToEnd = Vector3.Distance(testPoint, arcEnd);
                if (distanceToBegin < distanceToEnd)
                {
                    return arcStart;
                }
                else
                {
                    return arcEnd;
                }
            }
        }

        /// <summary>
        /// Compares two Vector3s and returns true if they are similar.
        /// </summary>
        /// <param name="vec1"></param>
        /// <param name="vec2"></param>
        /// <returns></returns>
        public static bool Approximately(Vector3 vec1, Vector3 vec2, float epsilon = 0.0001f)
        {
            return Vector3.Distance(vec1.normalized, vec2.normalized) < epsilon;
        }

        /// <summary>
        /// Given a direction, will return the closest cardinal-direction vector.
        /// </summary>
        /// <param name="direction"></param>
        /// <param name="magnitude"></param>
        /// <returns></returns>
        public static Vector3 GetCardinalDirection(Vector3 direction, out float magnitude)
        {
            float absX = Mathf.Abs(direction.x);
            float absY = Mathf.Abs(direction.y);
            float absZ = Mathf.Abs(direction.z);

            float dirX = direction.x / absX;
            float dirY = direction.y / absY;
            float dirZ = direction.z / absZ;

            if (absX > absY && absX > absZ)
            {
                magnitude = dirX;
                return new Vector3(dirX, 0, 0);
            }
            else if (absY > absX && absY > absZ)
            {
                magnitude = dirY;
                return new Vector3(0, dirY, 0);
            }
            else if (absZ > absX && absZ > absY)
            {
                magnitude = dirZ;
                return new Vector3(0, 0, dirZ);
            }
            else
            {
                magnitude = dirX;
                return new Vector3(dirX, 0, 0);
            }
        }

        /// <summary>
        /// Returns a Vector3 whose elements are all absolute values.
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Vector3 VectorAbs(Vector3 vector)
        {
            return new Vector3(Mathf.Abs(vector.x), Mathf.Abs(vector.y), Mathf.Abs(vector.z));
        }

        /// <summary>
        /// Returns a Vector2 whose elements are all absolute values.
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Vector2 VectorAbs(Vector2 vector)
        {
            return new Vector2(Mathf.Abs(vector.x), Mathf.Abs(vector.y));
        }
        #endregion


        #region Intergrators
        /// <summary>
        /// Simple verlet integrator.
        /// Often used for stepping physics updates.
        /// </summary>
        /// <param name="dt">Dt.</param>
        public static void Verlet(float dt, ref Vector2 pos, ref Vector2 vel)
        {
            Vector2 delta = (vel * dt) + 0.5f * (Physics2D.gravity * dt * dt);
            pos += delta;
            vel += Physics2D.gravity * dt;
            //return new Vector2[] { pos, vel };
        }

        /// <summary>
        /// Simple vertlet integrator.
        /// Often used for stepping physics updates.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="gravity"></param>
        /// <param name="pos"></param>
        /// <param name="vel"></param>
        public static void Verlet(float dt, Vector3 gravity, ref Vector3 pos, ref Vector3 vel)
        {
            Vector3 delta = (vel * dt) + 0.5f * (gravity * dt * dt);
            pos += delta;
            vel += gravity * dt;
        }

        /// <summary>
        /// Simple verlet integrator.
        /// Often used for stepping physics updates.
        /// </summary>
        /// <param name="dt">Dt.</param>
        public static void Verlet(float dt, float gravity, ref float pos, ref float vel)
        {
            float delta = (vel * dt) + 0.5f * (gravity * dt * dt);
            pos += delta;
            vel += gravity * dt;
        }

        public struct State
        {
            public float x;          // position
            public float v;          // velocity

            public State(float x, float v)
            {
                this.x = x;
                this.v = v;
            }
        };

        public struct Derivative
        {
            public float dx;          // derivative of position: velocity
            public float dv;          // derivative of velocity: acceleration

            public Derivative(float x, float v)
            {
                dx = x;
                dv = v;
            }
        };

        public static Derivative Evaluate(State initial, float t, float dt, Derivative d)
        {
            State state = new State(initial.x + d.dx * dt, initial.v + d.dv * dt);

            Derivative output = new Derivative();
            output.dx = state.v;
            output.dv = Acceleration(state, t + dt);
            return output;
        }

        public static float Acceleration(State state, float t)
        {
            const float k = 10;
            const float b = 1;
            return -k * state.x - b * state.v;
        }

        /// <summary>
        /// Runge-Kutta four-term intergator.
        /// Often used for stepping physics updates.
        /// </summary>
        /// <param name="state"></param>
        /// <param name="t"></param>
        /// <param name="dt"></param>
        /// <returns></returns>
        public static State RK4Integrate(State state, float t, float dt)
        {
            Derivative a = Evaluate(state, t, 0.0f, new Derivative());
            Derivative b = Evaluate(state, t, dt * 0.5f, a);
            Derivative c = Evaluate(state, t, dt * 0.5f, b);
            Derivative d = Evaluate(state, t, dt, c);

            float dxdt = 1.0f / 6.0f * (a.dx + 2.0f * (b.dx + c.dx) + d.dx);
            float dvdt = 1.0f / 6.0f * (a.dv + 2.0f * (b.dv + c.dv) + d.dv);

            return new State(state.x + dxdt * dt, state.v + dvdt * dt);
        }

        /// <summary>
        /// Runge-Kutta four-term intergator.
        /// Often used for stepping physics updates.
        /// </summary>
        /// <param name="state"></param>
        /// <param name="t"></param>
        /// <param name="dt"></param>
		public static void RK4Integrate(ref State state, float t, float dt)
        {
            Derivative a = Evaluate(state, t, 0.0f, new Derivative());
            Derivative b = Evaluate(state, t, dt * 0.5f, a);
            Derivative c = Evaluate(state, t, dt * 0.5f, b);
            Derivative d = Evaluate(state, t, dt, c);

            float dxdt = 1.0f / 6.0f * (a.dx + 2.0f * (b.dx + c.dx) + d.dx);
            float dvdt = 1.0f / 6.0f * (a.dv + 2.0f * (b.dv + c.dv) + d.dv);

            state.x = state.x + dxdt * dt;
            state.v = state.v + dvdt * dt;
        }
        #endregion


        #region Positional
        /// <summary>
        /// 
        /// </summary>
        /// <param name="targetDir"></param>
        /// <param name="angle"></param>
        /// <returns></returns>
        public static Vector3 GetPointOnUnitSphereCap(Quaternion targetDir, float angle)
        {
            var angledInRad = UnityEngine.Random.Range(0.0f, angle) * Mathf.Deg2Rad;
            var PointOnCircle = (UnityEngine.Random.insideUnitCircle.normalized) * Mathf.Sin(angledInRad);
            var vec = new Vector3(PointOnCircle.x, PointOnCircle.y, Mathf.Cos(angledInRad));
            return targetDir * vec;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="targetDir"></param>
        /// <param name="angle"></param>
        /// <returns></returns>
        public static Vector3 GetPointOnUnitSphereCap(Vector3 targetDir, float angle)
        {
            return GetPointOnUnitSphereCap(Quaternion.LookRotation(targetDir), angle);
        }

        /// <summary>
        /// Returns a point in world space that would move a transform's parent so that the child centered on the given worldPoint.
        /// </summary>
        /// <param name="anchor"></param>
        /// <param name="offset"></param>
        /// <returns></returns>
        public static Vector3 GetRelativeWorldPosition(Vector3 worldPoint, Transform child)
        {
            Transform parent = child.transform;
            Assert.IsNotNull(parent);
            return parent.position - child.localPosition;
        }

        /// <summary>
        /// Rotates and positions the parent so that the child is facing against the anchor.
        /// </summary>
        void MakeChildFacingAnchor(Transform parent, Transform child, Transform anchor)
        {
            parent.rotation = Quaternion.LookRotation(-anchor.forward, anchor.up) * Quaternion.Inverse(Quaternion.Inverse(parent.rotation) * child.rotation);
            parent.position = anchor.position - (child.position - parent.position);
        }

        static List<Vector3> TempVecList = new List<Vector3>(10);
        static HashSet<Vector3> TempVecSet = new HashSet<Vector3>();
        /// <summary>
        /// Returns a volatile list of collider centers. This list is a shared value and should be considered temporary.
        /// If it must be cached, be sure to make a copy of the contents.
        /// </summary>
        /// <param name="cols"></param>
        /// <param name="excludeSelf"></param>
        /// <returns></returns>
        public static List<Vector3> GetColliderCenters(Collider[] cols)
        {
            TempVecSet.Clear();
            TempVecList.Clear();

            for (int i = 0; i < cols.Length; i++)
                TempVecSet.Add(cols[i].bounds.center);

            foreach (var pos in TempVecSet)
                TempVecList.Add(pos);

            return TempVecList;
        }

        /// <summary>
        /// Returns a volatile list of collider centers. This list is a shared value and should be considered temporary.
        /// If it must be cached, be sure to make a copy of the contents.
        /// </summary>
        /// <param name="cols"></param>
        /// <param name="excludeSelf"></param>
        /// <param name="len"></param>
        /// <returns></returns>
        public static List<Vector3> GetColliderCenters(Collider[] cols, int len)
        {
            TempVecSet.Clear();
            TempVecList.Clear();

            for (int i = 0; i < len; i++)
                TempVecSet.Add(cols[i].bounds.center);

            foreach (var pos in TempVecSet)
                TempVecList.Add(pos);

            return TempVecList;
        }

        static List<Transform> TempTransList = new List<Transform>(10);
        static HashSet<Transform> TempTransSet = new HashSet<Transform>();
        /// <summary>
        /// 
        /// </summary>
        /// <param name="cols"></param>
        /// <param name="excludeSelf"></param>
        /// <returns>A volitile list of Transforms representing the centers of each collider.</returns>
        public static List<Transform> GetColliderTransforms(Collider[] cols)
        {
            TempTransSet.Clear();
            TempTransList.Clear();

            for (int i = 0; i < cols.Length; i++)
                TempTransSet.Add(cols[i].transform);

            foreach (var pos in TempTransSet)
                TempTransList.Add(pos);

            return TempTransList;
        }

        /// <summary>
        /// Similar to GetColliderTransforms but it takes a length representing the number of valid entries.
        /// </summary>
        /// <param name="cols"></param>
        /// <returns></returns>
        public static List<Transform> GetColliderTransforms(Collider[] cols, int len)
        {
            TempTransSet.Clear();
            TempTransList.Clear();

            for (int i = 0; i < len; i++)
                TempTransSet.Add(cols[i].transform);

            foreach (var pos in TempTransSet)
                TempTransList.Add(pos);

            return TempTransList;
        }

        /// <summary>
        /// Returns the transform in the list that is closest to the given point.
        /// </summary>
        /// <param name="point"></param>
        public static Transform GetClosest(List<Transform> trans, Vector3 point)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);

            int bestIndex = 0;
            float bestMatchMag = float.MaxValue;
            for (int i = 0; i < trans.Count; i++)
            {
                float mag = (trans[i].transform.position - point).sqrMagnitude;
                if (mag <= bestMatchMag)
                {
                    bestMatchMag = mag;
                    bestIndex = i;
                }
            }

            return trans[bestIndex];
        }
        
        /// <summary>
        /// Returns the point in the list that is closest to the given point.
        /// </summary>
        /// <param name="point"></param>
        public static Vector3 GetClosest(List<Vector3> positions, Vector3 point)
        {
            Assert.IsNotNull(positions);
            Assert.IsTrue(positions.Count > 0);

            int bestIndex = 0;
            float bestMatchMag = float.MaxValue;
            for (int i = 0; i < positions.Count; i++)
            {
                float mag = (positions[i] - point).sqrMagnitude;
                if (mag <= bestMatchMag)
                {
                    bestMatchMag = mag;
                    bestIndex = i;
                }
            }

            return positions[bestIndex];
        }

        /// <summary>
        /// Returns the index of the transform in the list that is closest to the given point.
        /// </summary>
        /// <param name="point"></param>
        public static int GetIndexOfClosest(List<Transform> trans, Vector3 point)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);

            int bestIndex = 0;
            float bestMatchMag = float.MaxValue;
            for (int i = 0; i < trans.Count; i++)
            {
                float mag = (trans[i].transform.position - point).sqrMagnitude;
                if (mag <= bestMatchMag)
                {
                    bestMatchMag = mag;
                    bestIndex = i;
                }
            }

            return bestIndex;
        }



        /// <summary>
        /// Returns the index of the point in the list that is closest to the given point.
        /// </summary>
        /// <param name="point"></param>
        public static int GetIndexOfClosest(List<Vector3> positions, Vector3 point)
        {
            Assert.IsNotNull(positions);
            Assert.IsTrue(positions.Count > 0);

            int bestIndex = 0;
            float bestMatchMag = float.MaxValue;
            for (int i = 0; i < positions.Count; i++)
            {
                float mag = (positions[i] - point).sqrMagnitude;
                if (mag <= bestMatchMag)
                {
                    bestMatchMag = mag;
                    bestIndex = i;
                }
            }

            return bestIndex;
        }

        /// <summary>
        /// Returns the index of the RaycastHit in the list whose contact point is closest to the given point.
        /// </summary>
        /// <param name="center"></param>
        public static int GetIndexOfClosest(int hitCount, RaycastHit[] hits, Vector3 center)
        {
            Assert.IsNotNull(hits);
            Assert.IsTrue(hitCount > 0);

            int bestIndex = 0;
            float bestMatchMag = float.MaxValue;
            for (int i = 0; i < hitCount; i++)
            {
                float mag = (hits[i].point - center).sqrMagnitude;
                if (mag <= bestMatchMag)
                {
                    bestMatchMag = mag;
                    bestIndex = i;
                }
            }

            return bestIndex;
        }

        /// <summary>
        /// Returns the transform in the list that is furthest from the given point.
        /// </summary>
        /// <param name="point"></param>
        public static Transform GetFurthest(List<Transform> trans, Vector3 point)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);

            int bestIndex = 0;
            float bestMatchMag =  -float.MaxValue;
            for (int i = 0; i < trans.Count; i++)
            {
                float mag = (trans[i].transform.position - point).sqrMagnitude;
                if (mag >= bestMatchMag)
                {
                    bestMatchMag = mag;
                    bestIndex = i;
                }
            }

            return trans[bestIndex];
        }

        /// <summary>
        /// Returns the point in the list that is furthest from the given point.
        /// </summary>
        /// <param name="point"></param>
        public static Vector3 GetFurthest(List<Vector3> trans, Vector3 point)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);

            int bestIndex = 0;
            float bestMatchMag = -float.MaxValue;
            for (int i = 0; i < trans.Count; i++)
            {
                float mag = (trans[i] - point).sqrMagnitude;
                if (mag >= bestMatchMag)
                {
                    bestMatchMag = mag;
                    bestIndex = i;
                }
            }

            return trans[bestIndex];
        }

        /// <summary>
        /// Returns the index of the RaycastHit in the list whose contact point is furthest from the given point.
        /// </summary>
        /// <param name="center"></param>
        public static int GetIndexOfFurthest(int hitCount, RaycastHit[] hits, Vector3 center)
        {
            Assert.IsNotNull(hits);
            Assert.IsTrue(hitCount > 0);

            int bestIndex = 0;
            float bestMatchMag = -float.MaxValue;
            for (int i = 0; i < hitCount; i++)
            {
                float mag = (hits[i].point - center).sqrMagnitude;
                if (mag >= bestMatchMag)
                {
                    bestMatchMag = mag;
                    bestIndex = i;
                }
            }

            return bestIndex;
        }

        /// <summary>
        /// Sorts a list of transforms closest-to-furthest.
        /// </summary>
        /// <param name="trans"></param>
        public static void SortClosestFirst(Vector3 p, List<Transform> trans)
        {
            trans.Sort(CompareTransformDistance.CompareDist(p));
            //trans.Sort((a, b) => { return (a.transform.position - p).sqrMagnitude.CompareTo((b.transform.position - p).sqrMagnitude); });
        }
        sealed class CompareTransformDistance : IComparer<Transform>
        {
            static CompareTransformDistance Shared = new CompareTransformDistance();
            public static IComparer<Transform> CompareDist(Vector3 p)
            {
                Shared.p = p;
                return Shared;
            }
            Vector3 p;
            public int Compare(Transform a, Transform b)
            {
                return (a.transform.position - p).sqrMagnitude.CompareTo((b.transform.position - p).sqrMagnitude);
            }
        }

        /// <summary>
        /// Sorts a list of transforms closest-to-furthest.
        /// </summary>
        /// <param name="trans"></param>
        public static void SortClosestFirst(Vector3 p, List<Vector3> positions)
        {
            positions.Sort(ComparePositionDistance.CompareDist(p));
        }
        sealed class ComparePositionDistance : IComparer<Vector3>
        {
            static ComparePositionDistance Shared = new ComparePositionDistance();
            public static IComparer<Vector3> CompareDist(Vector3 p)
            {
                Shared.p = p;
                return Shared;
            }
            Vector3 p;
            public int Compare(Vector3 a, Vector3 b)
            {
                return (a - p).sqrMagnitude.CompareTo((b - p).sqrMagnitude);
            }
        }

        static float xMin, yMin, zMin, xMax, yMax, zMax;
        /// <summary>
        /// Returns the geometric center of a collection of Transforms.
        /// </summary>
        /// <param name="trans"></param>
        /// <returns></returns>
        public static Vector3 GetCentroid(List<Transform> trans)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);

            var min = trans[0].position;
            var max = min;
            xMin = min.x;
            yMin = min.y;
            zMin = min.z;

            xMax = max.x;
            yMax = max.y;
            zMax = max.z;

            for (int i = 1; i < trans.Count; i++)
            {
                Vector3 pos = trans[i].position;
                if (pos.x < xMin) xMin = pos.x;
                if (pos.y < yMin) yMin = pos.y;
                if (pos.z < zMin) zMin = pos.z;

                if (pos.x > xMax) xMax = pos.x;
                if (pos.y > yMax) yMax = pos.y;
                if (pos.z > zMax) zMax = pos.z;
            }

            min = new Vector3(xMin, yMin, zMin);
            max = new Vector3(xMax, yMax, zMax);
            return min + ((max - min) * 0.5f);
        }

        static List<Vector4> SharedPoints = new List<Vector4>(20);
        /// <summary>
        /// Returns the geometric center of a collection of Transforms.
        /// </summary>
        /// <param name="trans"></param>
        /// <returns></returns>
        public static Vector3 GetCentroid(List<Transform> trans, List<float> weights)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);
            Assert.IsTrue(trans.Count == weights.Count);
            
            int count = trans.Count;
            SharedPoints.Clear();
            for(int i = 0; i < count; i++)
            {
                var pos = trans[i].position;
                SharedPoints.Add(new Vector4(pos.x, pos.y, pos.z, weights[i]));
            }

            return GetWeightedAABB(SharedPoints);
        }

        /// <summary>
        /// Returns the geometric center of a collection of Transforms.
        /// </summary>
        /// <param name="trans"></param>
        /// <returns></returns>
        public static Vector3 GetLimitedCentroid(List<Transform> trans, List<float> weights)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);
            Assert.IsTrue(trans.Count == weights.Count);
            
            int count = trans.Count;
            SharedPoints.Clear();
            for (int i = 0; i < count; i++)
            {
                var pos = trans[i].position;
                SharedPoints.Add(new Vector4(pos.x, pos.y, pos.z, weights[i]));
            }
            return GetWeightLimitedAABB(SharedPoints);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="points"></param>
        /// <param name="weightedCenter"></param>
        public static Vector3 GetWeightedAABB(Vector4[] points)
        {
            Vector3 min = new Vector3(float.MaxValue, float.MaxValue, float.MaxValue);
            Vector3 max = new Vector3(float.MinValue, float.MinValue, float.MinValue);
            Vector3 minW = new Vector3(float.MinValue, float.MinValue, float.MinValue);
            Vector3 maxW = new Vector3(float.MinValue, float.MinValue, float.MinValue);

            for (int i = 0; i < points.Length; i++)
            {
                TrySetExtent((Vector3)points[i] - min, points[i], ref min, ref minW);
                TrySetExtent(max - (Vector3)points[i], points[i], ref max, ref maxW);
            }

            var bounds = new Bounds();
            bounds.SetMinMax(min, max);
            Vector3 center_b = bounds.center;
            Vector3 sumW = minW + maxW;

            var weightedCenter = new Vector3()
            {
                x = Mathf.Abs(sumW.x) < float.Epsilon ? center_b.x : Mathf.LerpUnclamped(min.x, max.x, maxW.x / sumW.x),
                y = Mathf.Abs(sumW.y) < float.Epsilon ? center_b.y : Mathf.LerpUnclamped(min.y, max.y, maxW.y / sumW.y),
                z = Mathf.Abs(sumW.z) < float.Epsilon ? center_b.z : Mathf.LerpUnclamped(min.z, max.z, maxW.z / sumW.z)
            };

            return weightedCenter;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="points"></param>
        /// <param name="weightedCenter"></param>
        public static Vector3 GetWeightedAABB(List<Vector4> points)
        {
            Vector3 min = new Vector3(float.MaxValue, float.MaxValue, float.MaxValue);
            Vector3 max = new Vector3(float.MinValue, float.MinValue, float.MinValue);
            Vector3 minW = new Vector3(float.MinValue, float.MinValue, float.MinValue);
            Vector3 maxW = new Vector3(float.MinValue, float.MinValue, float.MinValue);

            int count = points.Count;
            for (int i = 0; i < count; i++)
            {
                TrySetExtent((Vector3)points[i] - min, points[i], ref min, ref minW);
                TrySetExtent(max - (Vector3)points[i], points[i], ref max, ref maxW);
            }

            var bounds = new Bounds();
            bounds.SetMinMax(min, max);
            Vector3 center_b = bounds.center;
            Vector3 sumW = minW + maxW;

            var weightedCenter = new Vector3()
            {
                x = Mathf.Abs(sumW.x) < float.Epsilon ? center_b.x : Mathf.LerpUnclamped(min.x, max.x, maxW.x / sumW.x),
                y = Mathf.Abs(sumW.y) < float.Epsilon ? center_b.y : Mathf.LerpUnclamped(min.y, max.y, maxW.y / sumW.y),
                z = Mathf.Abs(sumW.z) < float.Epsilon ? center_b.z : Mathf.LerpUnclamped(min.z, max.z, maxW.z / sumW.z)
            };

            return weightedCenter;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="points"></param>
        /// <param name="weightedCenter"></param>
        public static Vector3 GetWeightLimitedAABB(Vector4[] points)
        {
            var trueCentroid = GetCentroid4(points);

            //now we need to calculate each point's 'pull' away from the true centroid
            float maxPull = Mathf.Epsilon;
            int len = points.Length;

            for (int i = 0; i < len; i++)
                maxPull += points[i].w;

            float diff;
            Vector4 p;
            Vector4 composite = trueCentroid;
            for(int i = 0; i < len; i++)
            {
                p = points[i];
                diff = p.w / maxPull;
                Vector4 vec = p - trueCentroid;
                composite += (vec * diff);
            }

            return composite;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="points"></param>
        /// <param name="weightedCenter"></param>
        public static Vector3 GetWeightLimitedAABB(List<Vector4> points)
        {
            var trueCentroid = GetCentroid4(points);

            //now we need to calculate each point's 'pull' away from the true centroid
            float maxPull = Mathf.Epsilon;
            int len = points.Count;

            for (int i = 0; i < len; i++)
                maxPull += points[i].w;

            float diff;
            Vector4 p;
            Vector4 composite = trueCentroid;
            for (int i = 0; i < len; i++)
            {
                p = points[i];
                diff = p.w / maxPull;
                Vector4 vec = p - trueCentroid;
                composite += (vec * diff);
            }

            return composite;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="test"></param>
        /// <param name="point"></param>
        /// <param name="limit"></param>
        /// <param name="limitWeight"></param>
        private static void TrySetExtent(Vector3 test, Vector4 point, ref Vector3 limit, ref Vector3 limitWeight)
        {
            if (test.x < 0 || (test.x < float.Epsilon && point.w > limitWeight.x))
            {
                limit.x = point.x;
                limitWeight.x = point.w;
            }

            if (test.y < 0 || (test.y < float.Epsilon && point.w > limitWeight.y))
            {
                limit.y = point.y;
                limitWeight.y = point.w;
            }

            if (test.z < 0 || (test.z < float.Epsilon && point.w > limitWeight.z))
            {
                limit.z = point.z;
                limitWeight.z = point.w;
            }
        }

        static List<Vector3> TempVec3s = new List<Vector3>(5);
        /// <summary>
        /// Returns the geometric center of a collection of Transforms and Rigidbodies.
        /// </summary>
        /// <param name="trans"></param>
        /// <param name="bodies"></param>
        /// <returns></returns>
        public static Vector3 GetCentroid(List<Transform> trans, List<Rigidbody> bodies)
        {
            //var TempVec3s = SharedArrayFactory.RequestTempList<Vector3>();
            TempVec3s.Clear();
            for (int i = 0; i < trans.Count; i++)
                TempVec3s.Add(trans[i].position);
            for (int i = 0; i < bodies.Count; i++)
                TempVec3s.Add(bodies[i].position);
            return MathUtils.GetCentroid(TempVec3s);
        }

        /// <summary>
        /// Returns the geometric center of a collection of Vector3s.
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static Vector3 GetCentroid(List<Vector3> points)
        {
            Assert.IsNotNull(points);
            Assert.IsTrue(points.Count > 0);

            var min = points[0];
            var max = min;
            xMin = min.x;
            yMin = min.y;
            zMin = min.z;

            xMax = max.x;
            yMax = max.y;
            zMax = max.z;

            int len = points.Count;
            for (int i = 1; i < len; i++)
            {
                Vector3 pos = points[i];
                if (pos.x < xMin) xMin = pos.x;
                if (pos.y < yMin) yMin = pos.y;
                if (pos.z < zMin) zMin = pos.z;

                if (pos.x > xMax) xMax = pos.x;
                if (pos.y > yMax) yMax = pos.y;
                if (pos.z > zMax) zMax = pos.z;
            }

            min = new Vector3(xMin, yMin, zMin);
            max = new Vector3(xMax, yMax, zMax);
            return min + ((max - min) * 0.5f);
        }

        /// <summary>
        /// Returns the geometric center of a collection of Vector3s.
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static Vector3 GetCentroid(List<Vector4> points)
        {
            Assert.IsNotNull(points);
            Assert.IsTrue(points.Count > 0);

            var min = points[0];
            var max = min;
            xMin = min.x;
            yMin = min.y;
            zMin = min.z;

            xMax = max.x;
            yMax = max.y;
            zMax = max.z;

            int len = points.Count;
            for (int i = 1; i < len; i++)
            {
                Vector3 pos = points[i];
                if (pos.x < xMin) xMin = pos.x;
                if (pos.y < yMin) yMin = pos.y;
                if (pos.z < zMin) zMin = pos.z;

                if (pos.x > xMax) xMax = pos.x;
                if (pos.y > yMax) yMax = pos.y;
                if (pos.z > zMax) zMax = pos.z;
            }

            min = new Vector3(xMin, yMin, zMin);
            max = new Vector3(xMax, yMax, zMax);
            return min + ((max - min) * 0.5f);
        }

        /// <summary>
        /// Returns the geometric center of a collection of Vector3s.
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static Vector3 GetCentroid(Vector3[] points)
        {
            Assert.IsNotNull(points);
            Assert.IsTrue(points.Length > 0);

            var min = points[0];
            var max = min;
            xMin = min.x;
            yMin = min.y;
            zMin = min.z;

            xMax = max.x;
            yMax = max.y;
            zMax = max.z;

            for (int i = 1; i < points.Length; i++)
            {
                Vector3 pos = points[i];
                if (pos.x < xMin) xMin = pos.x;
                if (pos.y < yMin) yMin = pos.y;
                if (pos.z < zMin) zMin = pos.z;

                if (pos.x > xMax) xMax = pos.x;
                if (pos.y > yMax) yMax = pos.y;
                if (pos.z > zMax) zMax = pos.z;
            }

            min = new Vector3(xMin, yMin, zMin);
            max = new Vector3(xMax, yMax, zMax);
            return min + ((max - min) * 0.5f);
        }

        /// <summary>
        /// Returns the geometric center of a collection of Vector3s.
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static Vector3 GetCentroid(Vector4[]points)
        {
            Assert.IsNotNull(points);
            Assert.IsTrue(points.Length > 0);

            var min = points[0];
            var max = min;
            xMin = min.x;
            yMin = min.y;
            zMin = min.z;

            xMax = max.x;
            yMax = max.y;
            zMax = max.z;

            for (int i = 1; i < points.Length; i++)
            {
                Vector3 pos = points[i];
                if (pos.x < xMin) xMin = pos.x;
                if (pos.y < yMin) yMin = pos.y;
                if (pos.z < zMin) zMin = pos.z;

                if (pos.x > xMax) xMax = pos.x;
                if (pos.y > yMax) yMax = pos.y;
                if (pos.z > zMax) zMax = pos.z;
            }

            min = new Vector3(xMin, yMin, zMin);
            max = new Vector3(xMax, yMax, zMax);
            return min + ((max - min) * 0.5f);
        }

        /// <summary>
        /// Returns the geometric center of a collection of Vector3s.
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static Vector4 GetCentroid4(Vector4[] points)
        {
            Assert.IsNotNull(points);
            Assert.IsTrue(points.Length > 0);

            var min = points[0];
            var max = min;
            xMin = min.x;
            yMin = min.y;
            zMin = min.z;

            xMax = max.x;
            yMax = max.y;
            zMax = max.z;

            for (int i = 1; i < points.Length; i++)
            {
                Vector3 pos = points[i];
                if (pos.x < xMin) xMin = pos.x;
                if (pos.y < yMin) yMin = pos.y;
                if (pos.z < zMin) zMin = pos.z;

                if (pos.x > xMax) xMax = pos.x;
                if (pos.y > yMax) yMax = pos.y;
                if (pos.z > zMax) zMax = pos.z;
            }

            min = new Vector3(xMin, yMin, zMin);
            max = new Vector3(xMax, yMax, zMax);
            return min + ((max - min) * 0.5f);
        }

        /// <summary>
        /// Returns the geometric center of a collection of Vector3s.
        /// </summary>
        /// <param name="points"></param>
        /// <returns></returns>
        public static Vector4 GetCentroid4(List<Vector4> points)
        {
            Assert.IsNotNull(points);
            Assert.IsTrue(points.Count > 0);

            var min = points[0];
            var max = min;
            xMin = min.x;
            yMin = min.y;
            zMin = min.z;

            xMax = max.x;
            yMax = max.y;
            zMax = max.z;

            int len = points.Count;
            for (int i = 1; i < len; i++)
            {
                Vector3 pos = points[i];
                if (pos.x < xMin) xMin = pos.x;
                if (pos.y < yMin) yMin = pos.y;
                if (pos.z < zMin) zMin = pos.z;

                if (pos.x > xMax) xMax = pos.x;
                if (pos.y > yMax) yMax = pos.y;
                if (pos.z > zMax) zMax = pos.z;
            }

            min = new Vector3(xMin, yMin, zMin);
            max = new Vector3(xMax, yMax, zMax);
            return min + ((max - min) * 0.5f);
        }

        /// <summary>
        /// Returns the averaged position of a list of Transforms.
        /// </summary>
        /// <param name="group"></param>
        /// <returns></returns>
        public static Vector3 GetAveragePosition(List<Transform> trans)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);

            var avg = trans[0].position;
            for (int i = 1; i < trans.Count; i++)
                avg += trans[i].position;

            return avg / trans.Count;
        }

        /// <summary>
        /// Returns the averaged position of a list of Vector3s.
        /// </summary>
        /// <param name="group"></param>
        /// <returns></returns>
        public static Vector3 GetAveragePosition(List<Vector3> trans)
        {
            Assert.IsNotNull(trans);
            Assert.IsTrue(trans.Count > 0);

            var avg = trans[0];
            for (int i = 1; i < trans.Count; i++)
                avg += trans[i];

            return avg / trans.Count;
        }
        #endregion

        #region Motion

        public enum MovementTypes
        {
            Othagonal,
            Oblique,
        }

        /// <summary>
        /// Applies rotation to a movement vector so that the z-axis motion always appears vertical to a facing-target.
        /// </summary>
        /// <param name="vel"></param>
        /// <param name="faceTarget"></param>
        /// <returns></returns>
        public static Vector3 TransformByFacingSpace(Vector3 vel, Transform faceTarget, MovementTypes type)
        {
            //translate to camera relative motion
            if (vel.sqrMagnitude > 0.01f)
            {
                if (type == MovementTypes.Othagonal)
                    return (Quaternion.AngleAxis(faceTarget.eulerAngles.y, Vector3.up) * vel);
                else if (type == MovementTypes.Oblique)
                {
                    //this is trickier - we need rotate vertical motion and avoid rotating horizontal
                    //motion but ensure we /maintain the same total magnitude.
                    var hor = new Vector3(vel.x, 0);
                    var vert = new Vector3(0, 0, vel.z);
                    var comp = Quaternion.AngleAxis(faceTarget.eulerAngles.y, Vector3.up) * vert;
                    return (comp + hor).normalized * vel.magnitude;
                }
            }

            return vel;
        }

        /// <summary>
        /// Coroutine that can be used to move an object along a parabola-based jump arc over time.
        /// </summary>
        /// <param name="agent"></param>
        /// <param name="startPos"></param>
        /// <param name="endPos"></param>
        /// <param name="jumpHeight"></param>
        /// <param name="duration"></param>
        /// <returns></returns>
        public static IEnumerator MoveAlongParabola(Transform agent, Vector3 startPos, Vector3 endPos, float jumpHeight, float duration, Action callback = null)
        {
            float normalizedTime = 0.0f;
            while (normalizedTime < 1.0f)
            {
                float yOffset = jumpHeight * 4.0f * (normalizedTime - normalizedTime * normalizedTime);
                agent.position = Vector3.Lerp(startPos, endPos, normalizedTime) + yOffset * Vector3.up;
                normalizedTime += Time.deltaTime / duration;
                yield return null;
            }

            callback?.Invoke();
        }

        /// <summary>
        /// Coroutine that can be used to move an object along a curve-based jump arc over time.
        /// </summary>
        /// <param name="agent"></param>
        /// <param name="startPos"></param>
        /// <param name="endPos"></param>
        /// <param name="curve"></param>
        /// <param name="jumpHeight"></param>
        /// <param name="duration"></param>
        /// <returns></returns>
        public static IEnumerator MoveAlongCurve(Transform agent, Vector3 startPos, Vector3 endPos, AnimationCurve curve, float duration, Action callback = null)
        {
            float normalizedTime = 0.0f;
            while (normalizedTime < 1.0f)
            {
                float yOffset = curve.Evaluate(normalizedTime);
                agent.transform.position = Vector3.Lerp(startPos, endPos, normalizedTime) + yOffset * Vector3.up;
                normalizedTime += Time.deltaTime / duration;
                yield return null;
            }

            if (callback != null)
                callback();
        }
        #endregion


        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static int NearestSuperiorPowerOf2(int n)
        {
            return (int)Mathf.Pow(2, Mathf.Ceil(Mathf.Log(n) / Mathf.Log(2)));
        }

        /// <summary>
        /// Performs a fisher-yates shuffle on an array.
        /// The source array is left untouched. The out array
        /// cannot be the same as the in array.
        /// Both arrays must have the same size.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="listIn"></param>
        /// <param name="listOut"></param>
        /// <param name="rng">An optional RNG that can be used to provided seeded output. If null, UnityEngine.Random will be used instead.</param>
        public static void Shuffle<T>(T[] listIn, ref T[] listOut, System.Random rng = null)
        {
            Assert.IsNotNull(listIn);
            Assert.IsNotNull(listOut);
            Assert.AreEqual(listOut.Length, listIn.Length);
            Assert.AreNotEqual(listIn, listOut);

            int count = listIn.Length;
            if (count == 1)
            {
                listOut[0] = listIn[0];
                return;
            }

            for (var i = count - 1; i >= 0; i--)
            {
                var r = rng != null ? (float)rng.NextDouble() : UnityEngine.Random.value;
                var j = (int)Mathf.Floor(r * (i+1));
                listOut[i] = listIn[j];
                listOut[j] = listIn[i];
            }
        }

        /// <summary>
        /// Performs a fisher-yates shuffle on an array.
        /// The source array is left untouched. The out array
        /// cannot be the same as the in array.
        /// Both arrays must have the same size.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="listIn"></param>
        /// <param name="listOut"></param>
        /// <param name="rng">An optional RNG that can be used to provided seeded output. If null, UnityEngine.Random will be used instead.</param>
        public static void Shuffle<T>(List<T> listIn, ref List<T> listOut, System.Random rng = null)
        {
            Assert.IsNotNull(listIn);
            Assert.IsNotNull(listOut);
            Assert.AreEqual(listOut.Count, listIn.Count);
            Assert.AreNotEqual(listIn, listOut);

            int count = listIn.Count;
            if (count == 1)
            {
                listOut[0] = listIn[0];
                return;
            }

            for (var i = count - 1; i >= 0; i--)
            {
                var r = rng != null ? (float)rng.NextDouble() : UnityEngine.Random.value;
                var j = (int)Mathf.Floor(r * (i+1));
                listOut[i] = listIn[j];
                listOut[j] = listIn[i];
            }
        }

        
    }


    /// <summary>
    /// Extension methods for enumerables.
    /// </summary>
    public static class Enumerable
    {
        /// <summary>
        /// Performs a fisher-yates shuffle on an array.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="listIn"></param>
        /// <param name="listOut"></param>
        /// <param name="rng">An optional RNG that can be used to provided seeded output. If null, UnityEngine.Random will be used instead.</param>
        public static T[] Shuffle<T>(this T[] listIn, System.Random rng = null)
        {
            var listOut = new T[listIn.Length];
            MathUtils.Shuffle(listIn, ref listOut, rng);
            return listOut;
        }

        /// <summary>
        /// Performs a fisher-yates shuffle on a List<>.
        /// </summary>
        /// <typeparam name="TSource"></typeparam>
        /// <param name="listIn"></param>
        /// <param name="rng">An optional RNG that can be used to provided seeded output. If null, UnityEngine.Random will be used instead.</param>
        /// <returns></returns>
        public static List<TSource> Shuffle<TSource>(this List<TSource> listIn, System.Random rng = null)
        {
            var listOut = new List<TSource>(listIn.Count);
            MathUtils.Shuffle(listIn, ref listOut, rng);
            return listOut;
        }

        /// <summary>
        /// Performs a fisher-yates shuffle on a Stack<>.
        /// </summary>
        /// <typeparam name="TSource"></typeparam>
        /// <param name="listIn"></param>
        /// <param name="rng">An optional RNG that can be used to provided seeded output. If null, UnityEngine.Random will be used instead.</param>
        /// <returns></returns>
        public static Stack<TSource> Shuffle<TSource>(this Stack<TSource> stackIn, System.Random rng = null)
        {
            var stackOut = new Stack<TSource>();
            if (stackIn.Count == 1)
            {
                stackOut.Push(stackIn.Pop());
                return stackOut;
            }

            var listIn = stackIn.ToArray();
            var listOut = new TSource[listIn.Length];
            MathUtils.Shuffle(listIn, ref listOut, rng);

            for(int i = 0; i < listIn.Length; i++)
                stackOut.Push(listOut[0]);
            return stackOut;
        }

        /// <summary>
        /// Performs a fisher-yates shuffle on a Queue<>.
        /// </summary>
        /// <typeparam name="TSource"></typeparam>
        /// <param name="listIn"></param>
        /// <param name="rng">An optional RNG that can be used to provided seeded output. If null, UnityEngine.Random will be used instead.</param>
        /// <returns></returns>
        public static Queue<TSource> Shuffle<TSource>(this Queue<TSource> stackIn, System.Random rng = null)
        {
            var queueOut = new Queue<TSource>();
            if (stackIn.Count == 1)
            {
                queueOut.Enqueue(stackIn.Dequeue());
                return queueOut;
            }

            var listIn = stackIn.ToArray();
            var listOut = new TSource[listIn.Length];
            MathUtils.Shuffle(listIn, ref listOut, rng);

            for (int i = 0; i < listIn.Length; i++)
                queueOut.Enqueue(listOut[0]);
            return queueOut;
        }
    }

    
    /// <summary>
    /// 
    /// </summary>
    public static class LayerUtils
    {
        /// <summary>
        /// Checks to see if a LayerMask contains a layer index.
        /// </summary>
        /// <param name="mask"></param>
        /// <param name="layerIndex">A GameObject layer index. Must range between 0 and 31.</param>
        /// <returns></returns>
        public static bool ContainsLayerIndex(this LayerMask mask, int layerIndex)
        {
            Assert.IsTrue(layerIndex >= 0);
            Assert.IsTrue(layerIndex < 32);
            return (mask.value & (1 << layerIndex)) != 0;
        }

        /// <summary>
        /// Combines a list of GameObject layer indicies into a LayerMask that contains all of them.
        /// 
        /// NOT YET TESTED!
        /// 
        /// </summary>
        /// <param name="layerIndex"></param>
        /// <returns></returns>
        [Obsolete("Not yet tested!")]
        public static LayerMask Combine(params int[] layerIndex)
        {
            int value = 0;
            for (int i = 0; i < layerIndex.Length; i++)
            {
                Assert.IsTrue(layerIndex[i] >= 0);
                Assert.IsTrue(layerIndex[i] < 32);
                value |= (1 << layerIndex[i]);
            }
            var mask = new LayerMask();
            mask.value = value;
            return mask;
        }
    }


    /// <summary>
    /// Type that represents a value and a weight associated with it.
    /// Note that for this to appear in the Unity editor, concrete types
    /// must be derived from it.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    [System.Serializable]
    public class WeightedValue<T>
    {
        public T Value;
        public float Weight = 1;
    }

    /// <summary>
    /// A float value with a weight.
    /// </summary>
    [System.Serializable]
    public class WeightedFloat : WeightedValue<float> { }

    /// <summary>
    /// A float value with a weight.
    /// </summary>
    [System.Serializable]
    public class WeightedInt : WeightedValue<int> { }

    /// <summary>
    /// A float value with a weight.
    /// </summary>
    [System.Serializable]
    public class WeightedString : WeightedValue<string> { }


    /// <summary>
    /// Provides additional functionality to collections of weighted values.
    /// </summary>
    public static class WeightedValueCollections
    {
        /// <summary>
        /// Selects a random element from the array of weighted values, taking each value's weight into consideration.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="array">The array to select from.</param>
        /// <param name="randomizer">An optional randomaizer. Only set if a specific seed is required.</param>
        /// <returns></returns>
        public static T SelectWeightedRandom<T>(this WeightedValue<T>[] array, System.Random randomizer = null)
        {
            return array.SelectWeightedRandom(array.Sum(w => w.Weight), randomizer);
        }

        /// <summary>
        /// Selects a random element from the array of weighted values, taking each value's weight into consideration.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="array">The array to select from.</param>
        /// <param name="weightTotal">The combine total of all weights in the collection.  If this is known ahead of time, it can be supplied to save time not having to calculate it internally.</param>
        /// <param name="randomizer">An optional randomaizer. Only set if a specific seed is required.</param>
        /// <returns></returns>
        public static T SelectWeightedRandom<T>(this WeightedValue<T>[] array, double weightTotal, System.Random randomizer = null)
        {
            var r = randomizer == null ? UnityEngine.Random.Range(0, (float)weightTotal) : (randomizer.NextDouble() * weightTotal);
            double tally = r;
            for (int i = 0; i < array.Length; i++)
            {
                tally -= array[i].Weight;
                if (tally <= 0) return array[i].Value;
            }
            throw new UnityException("Weights not calculated correctly. Total: " + weightTotal + "  Random:  " + r + "  Final Tally: " + tally);
        }

        /// <summary>
        /// Selects a random element from the array of weighted values, taking each value's weight into consideration.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="collection">The collection to select from.</param>
        /// <param name="randomizer">An optional randomaizer. Only set if a specific seed is required.</param>
        /// <returns></returns>
        public static T SelectWeightedRandom<T>(this IEnumerable<WeightedValue<T>> collection, System.Random randomizer = null)
        {
            return collection.SelectWeightedRandom(collection.Sum(w => w.Weight), randomizer);
        }

        /// <summary>
        /// Selects a random element from the collection of weighted values, taking each value's weight into consideration.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="collection">The collection to choose from.</param>
        /// <param name="weightTotal">The combine total of all weights in the collection. If this is known ahead of time, it can be supplied to save time not having to calculate it internally.</param>
        /// <param name="randomizer">An optional randomizer. Only set if a specific seed is required.</param>
        /// <returns></returns>
        public static T SelectWeightedRandom<T>(this IEnumerable<WeightedValue<T>> collection, double weightTotal, System.Random randomizer = null)
        {
            var r = randomizer == null ? UnityEngine.Random.Range(0, (float)weightTotal) : (randomizer.NextDouble() * weightTotal);
            double tally = r;
            foreach(var element in collection)
            {
                tally -= element.Weight;
                if (tally <= 0) return element.Value;
            }
            throw new UnityException("Weights not calculated correctly. Total: " + weightTotal + "  Random:  " + r + "  Final Tally: " + tally);
        }
    }
    

    /// <summary>
    /// A concrete implementation of WeightedValue<T> for GameObjects.
    /// </summary>
    [System.Serializable]
    public class WeightedGameObject : WeightedValue<GameObject> { }
    [Serializable]
    public class WeightedBehaviour : WeightedValue<Behaviour> { }
    [Serializable]
    public class WeightedCollider : WeightedValue<Collider> { }
    [Serializable]
    public class WeightedCollider2D : WeightedValue<Collider2D> { }
    [Serializable]
    public class WeightedRenderer : WeightedValue<Renderer> { }
    [Serializable]
    public class WeightedRigidbody: WeightedValue<Rigidbody> { }
    [System.Serializable]
    public class WeightedColor : WeightedValue<Color> { }


    /// <summary>
    /// Provides extension methods for several built-in unity classes and structs.
    /// </summary>
    public static class MathExtensions
    {
        /// <summary>
        /// Returns true if two bound volumes intersect.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="aRot"></param>
        /// <param name="b"></param>
        /// <param name="bRot"></param>
        /// <returns></returns>
        public static bool Intersects(this Bounds a, Quaternion rot, Bounds otherBounds, Quaternion otheRot)
        {
            throw new UnityException("Not yet implemented.");
        }
    }
    

    /// <summary>
    /// 
    /// </summary>
    public static class TransformExtensions
    {
        /// <summary>
        /// Attempts to set a transform's scale in global-space. Due to the nature of lossy scaling, if the 
        /// object is parented there may be cases where the results will not be exact (due to rotations).
        /// </summary>
        /// <param name="transform"></param>
        /// <param name="globalScale"></param>
        public static void SetGlobalScale(this Transform transform, Vector3 globalScale)
        {
            transform.localScale = Vector3.one;
            transform.localScale = new Vector3(globalScale.x / transform.lossyScale.x, globalScale.y / transform.lossyScale.y, globalScale.z / transform.lossyScale.z);
        }
    }
    

    /// <summary>
    /// Extension methods for Unity's Vector3 nd Vector2 methods.
    /// </summary>
    public static class VectorExtensions
    {
        public enum SwizzelMode
        {
            XYZ = 0, //effectively, no change in most cases
            XZY,

            YXZ,
            YZX,

            ZXY,
            ZYX,
        }

        public enum SwizzelMode3DTo2D
        {
            XY,
            YX,
            XZ,
            ZX,
            YZ,
            ZY,
        }
        
        public enum SwizzelMode2D
        {
            XY,
            YX,
        }
        
        public static Vector3 Swizzel(this Vector3 vec, SwizzelMode mode)
        {
            switch (mode)
            {
                case SwizzelMode.XYZ:
                    {
                        return vec;
                    }
                case SwizzelMode.XZY:
                    {
                        return new Vector3(vec.x, vec.z, vec.y);
                    }
                case SwizzelMode.YXZ:
                    {
                        return new Vector3(vec.y, vec.x, vec.z);
                    }
                case SwizzelMode.YZX:
                    {
                        return new Vector3(vec.y, vec.z, vec.x);
                    }
                case SwizzelMode.ZXY:
                    {
                        return new Vector3(vec.z, vec.x, vec.y);
                    }
                case SwizzelMode.ZYX:
                    {
                        return new Vector3(vec.z, vec.y, vec.x);
                    }
                default:
                    {
                        return vec;
                    }
            }
        }
        
        public static Vector2 Swizzel(this Vector3 vec, SwizzelMode3DTo2D mode)
        {
            switch (mode)
            {
                case SwizzelMode3DTo2D.XY:
                    {
                        return vec;
                    }
                case SwizzelMode3DTo2D.XZ:
                    {
                        return new Vector2(vec.x, vec.z);
                    }
                case SwizzelMode3DTo2D.YX:
                    {
                        return new Vector2(vec.y, vec.x);
                    }
                case SwizzelMode3DTo2D.YZ:
                    {
                        return new Vector2(vec.y, vec.z);
                    }
                case SwizzelMode3DTo2D.ZX:
                    {
                        return new Vector2(vec.z, vec.x);
                    }
                case SwizzelMode3DTo2D.ZY:
                    {
                        return new Vector2(vec.z, vec.y);
                    }
                default:
                    {
                        return vec;
                    }
            }
        }
        
        public static Vector2 Swizzel(this Vector2 vec, SwizzelMode2D mode)
        {
            switch (mode)
            {
                case SwizzelMode2D.XY:
                    {
                        return vec;
                    }
                case SwizzelMode2D.YX:
                    {
                        return new Vector2(vec.y, vec.x);
                    }
                default:
                    {
                        return vec;
                    }
            }

        }
    }



    /// <summary>
    /// 
    /// </summary>
    public static class DataUtils
    {

        /// <summary>
        /// Encodes a single integer into a Vector3.
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public static Vector3 EncodeFloat3(int index)
        {
            float x = (float)(index & 1023) / 1023.0f;
            float y = (float)((index >> 10) & 1023) / 1023.0f;
            float z = (float)((index >> 20) & 511) / 511.0f;
            return new Vector3(x, y, z);
        }

        /// <summary>
        /// Decodes a Vector3 into a previously encoded integer.
        /// </summary>
        /// <param name="float3"></param>
        /// <returns></returns>
        public static int DecodeFloat3(Vector3 float3)
        {
            int x = (int)(float3.x * 1023);
            int y = (int)(float3.y * 1023) << 10;
            int z = (int)(float3.z * 511) << 20;
            return x | y | z;
        }
    }
}
 