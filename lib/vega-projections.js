(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(require('vega-projection')) :
  typeof define === 'function' && define.amd ? define(['vega-projection'], factory) :
  (factory(global.vega));
}(this, (function (vegaProjection) { 'use strict';

  // Adds floating point numbers with twice the normal precision.
  // Reference: J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and
  // Fast Robust Geometric Predicates, Discrete & Computational Geometry 18(3)
  // 305–363 (1997).
  // Code adapted from GeographicLib by Charles F. F. Karney,
  // http://geographiclib.sourceforge.net/

  function adder() {
    return new Adder;
  }

  function Adder() {
    this.reset();
  }

  Adder.prototype = {
    constructor: Adder,
    reset: function() {
      this.s = // rounded value
      this.t = 0; // exact error
    },
    add: function(y) {
      add(temp, y, this.t);
      add(this, temp.s, this.s);
      if (this.s) this.t += temp.t;
      else this.s = temp.t;
    },
    valueOf: function() {
      return this.s;
    }
  };

  var temp = new Adder;

  function add(adder, a, b) {
    var x = adder.s = a + b,
        bv = x - a,
        av = x - bv;
    adder.t = (a - av) + (b - bv);
  }

  var epsilon = 1e-6;
  var epsilon2 = 1e-12;
  var pi = Math.PI;
  var halfPi = pi / 2;
  var quarterPi = pi / 4;
  var tau = pi * 2;

  var degrees = 180 / pi;
  var radians = pi / 180;

  var abs = Math.abs;
  var atan = Math.atan;
  var atan2 = Math.atan2;
  var cos = Math.cos;
  var sin = Math.sin;
  var sqrt = Math.sqrt;

  function acos(x) {
    return x > 1 ? 0 : x < -1 ? pi : Math.acos(x);
  }

  function asin(x) {
    return x > 1 ? halfPi : x < -1 ? -halfPi : Math.asin(x);
  }

  function haversin(x) {
    return (x = sin(x / 2)) * x;
  }

  function noop() {}

  function streamGeometry(geometry, stream) {
    if (geometry && streamGeometryType.hasOwnProperty(geometry.type)) {
      streamGeometryType[geometry.type](geometry, stream);
    }
  }

  var streamObjectType = {
    Feature: function(object, stream) {
      streamGeometry(object.geometry, stream);
    },
    FeatureCollection: function(object, stream) {
      var features = object.features, i = -1, n = features.length;
      while (++i < n) streamGeometry(features[i].geometry, stream);
    }
  };

  var streamGeometryType = {
    Sphere: function(object, stream) {
      stream.sphere();
    },
    Point: function(object, stream) {
      object = object.coordinates;
      stream.point(object[0], object[1], object[2]);
    },
    MultiPoint: function(object, stream) {
      var coordinates = object.coordinates, i = -1, n = coordinates.length;
      while (++i < n) object = coordinates[i], stream.point(object[0], object[1], object[2]);
    },
    LineString: function(object, stream) {
      streamLine(object.coordinates, stream, 0);
    },
    MultiLineString: function(object, stream) {
      var coordinates = object.coordinates, i = -1, n = coordinates.length;
      while (++i < n) streamLine(coordinates[i], stream, 0);
    },
    Polygon: function(object, stream) {
      streamPolygon(object.coordinates, stream);
    },
    MultiPolygon: function(object, stream) {
      var coordinates = object.coordinates, i = -1, n = coordinates.length;
      while (++i < n) streamPolygon(coordinates[i], stream);
    },
    GeometryCollection: function(object, stream) {
      var geometries = object.geometries, i = -1, n = geometries.length;
      while (++i < n) streamGeometry(geometries[i], stream);
    }
  };

  function streamLine(coordinates, stream, closed) {
    var i = -1, n = coordinates.length - closed, coordinate;
    stream.lineStart();
    while (++i < n) coordinate = coordinates[i], stream.point(coordinate[0], coordinate[1], coordinate[2]);
    stream.lineEnd();
  }

  function streamPolygon(coordinates, stream) {
    var i = -1, n = coordinates.length;
    stream.polygonStart();
    while (++i < n) streamLine(coordinates[i], stream, 1);
    stream.polygonEnd();
  }

  function geoStream(object, stream) {
    if (object && streamObjectType.hasOwnProperty(object.type)) {
      streamObjectType[object.type](object, stream);
    } else {
      streamGeometry(object, stream);
    }
  }

  var areaRingSum = adder();

  var areaSum = adder(),
      lambda00,
      phi00,
      lambda0,
      cosPhi0,
      sinPhi0;

  var areaStream = {
    point: noop,
    lineStart: noop,
    lineEnd: noop,
    polygonStart: function() {
      areaRingSum.reset();
      areaStream.lineStart = areaRingStart;
      areaStream.lineEnd = areaRingEnd;
    },
    polygonEnd: function() {
      var areaRing = +areaRingSum;
      areaSum.add(areaRing < 0 ? tau + areaRing : areaRing);
      this.lineStart = this.lineEnd = this.point = noop;
    },
    sphere: function() {
      areaSum.add(tau);
    }
  };

  function areaRingStart() {
    areaStream.point = areaPointFirst;
  }

  function areaRingEnd() {
    areaPoint(lambda00, phi00);
  }

  function areaPointFirst(lambda, phi) {
    areaStream.point = areaPoint;
    lambda00 = lambda, phi00 = phi;
    lambda *= radians, phi *= radians;
    lambda0 = lambda, cosPhi0 = cos(phi = phi / 2 + quarterPi), sinPhi0 = sin(phi);
  }

  function areaPoint(lambda, phi) {
    lambda *= radians, phi *= radians;
    phi = phi / 2 + quarterPi; // half the angular distance from south pole

    // Spherical excess E for a spherical triangle with vertices: south pole,
    // previous point, current point.  Uses a formula derived from Cagnoli’s
    // theorem.  See Todhunter, Spherical Trig. (1871), Sec. 103, Eq. (2).
    var dLambda = lambda - lambda0,
        sdLambda = dLambda >= 0 ? 1 : -1,
        adLambda = sdLambda * dLambda,
        cosPhi = cos(phi),
        sinPhi = sin(phi),
        k = sinPhi0 * sinPhi,
        u = cosPhi0 * cosPhi + k * cos(adLambda),
        v = k * sdLambda * sin(adLambda);
    areaRingSum.add(atan2(v, u));

    // Advance the previous points.
    lambda0 = lambda, cosPhi0 = cosPhi, sinPhi0 = sinPhi;
  }

  function spherical(cartesian) {
    return [atan2(cartesian[1], cartesian[0]), asin(cartesian[2])];
  }

  function cartesian(spherical) {
    var lambda = spherical[0], phi = spherical[1], cosPhi = cos(phi);
    return [cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi)];
  }

  function cartesianDot(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  function cartesianCross(a, b) {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  }

  // TODO return a
  function cartesianAddInPlace(a, b) {
    a[0] += b[0], a[1] += b[1], a[2] += b[2];
  }

  function cartesianScale(vector, k) {
    return [vector[0] * k, vector[1] * k, vector[2] * k];
  }

  // TODO return d
  function cartesianNormalizeInPlace(d) {
    var l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    d[0] /= l, d[1] /= l, d[2] /= l;
  }

  var lambda0$1, phi0, lambda1, phi1, // bounds
      lambda2, // previous lambda-coordinate
      lambda00$1, phi00$1, // first point
      p0, // previous 3D point
      deltaSum = adder(),
      ranges,
      range;

  var boundsStream = {
    point: boundsPoint,
    lineStart: boundsLineStart,
    lineEnd: boundsLineEnd,
    polygonStart: function() {
      boundsStream.point = boundsRingPoint;
      boundsStream.lineStart = boundsRingStart;
      boundsStream.lineEnd = boundsRingEnd;
      deltaSum.reset();
      areaStream.polygonStart();
    },
    polygonEnd: function() {
      areaStream.polygonEnd();
      boundsStream.point = boundsPoint;
      boundsStream.lineStart = boundsLineStart;
      boundsStream.lineEnd = boundsLineEnd;
      if (areaRingSum < 0) lambda0$1 = -(lambda1 = 180), phi0 = -(phi1 = 90);
      else if (deltaSum > epsilon) phi1 = 90;
      else if (deltaSum < -epsilon) phi0 = -90;
      range[0] = lambda0$1, range[1] = lambda1;
    }
  };

  function boundsPoint(lambda, phi) {
    ranges.push(range = [lambda0$1 = lambda, lambda1 = lambda]);
    if (phi < phi0) phi0 = phi;
    if (phi > phi1) phi1 = phi;
  }

  function linePoint(lambda, phi) {
    var p = cartesian([lambda * radians, phi * radians]);
    if (p0) {
      var normal = cartesianCross(p0, p),
          equatorial = [normal[1], -normal[0], 0],
          inflection = cartesianCross(equatorial, normal);
      cartesianNormalizeInPlace(inflection);
      inflection = spherical(inflection);
      var delta = lambda - lambda2,
          sign$$1 = delta > 0 ? 1 : -1,
          lambdai = inflection[0] * degrees * sign$$1,
          phii,
          antimeridian = abs(delta) > 180;
      if (antimeridian ^ (sign$$1 * lambda2 < lambdai && lambdai < sign$$1 * lambda)) {
        phii = inflection[1] * degrees;
        if (phii > phi1) phi1 = phii;
      } else if (lambdai = (lambdai + 360) % 360 - 180, antimeridian ^ (sign$$1 * lambda2 < lambdai && lambdai < sign$$1 * lambda)) {
        phii = -inflection[1] * degrees;
        if (phii < phi0) phi0 = phii;
      } else {
        if (phi < phi0) phi0 = phi;
        if (phi > phi1) phi1 = phi;
      }
      if (antimeridian) {
        if (lambda < lambda2) {
          if (angle(lambda0$1, lambda) > angle(lambda0$1, lambda1)) lambda1 = lambda;
        } else {
          if (angle(lambda, lambda1) > angle(lambda0$1, lambda1)) lambda0$1 = lambda;
        }
      } else {
        if (lambda1 >= lambda0$1) {
          if (lambda < lambda0$1) lambda0$1 = lambda;
          if (lambda > lambda1) lambda1 = lambda;
        } else {
          if (lambda > lambda2) {
            if (angle(lambda0$1, lambda) > angle(lambda0$1, lambda1)) lambda1 = lambda;
          } else {
            if (angle(lambda, lambda1) > angle(lambda0$1, lambda1)) lambda0$1 = lambda;
          }
        }
      }
    } else {
      ranges.push(range = [lambda0$1 = lambda, lambda1 = lambda]);
    }
    if (phi < phi0) phi0 = phi;
    if (phi > phi1) phi1 = phi;
    p0 = p, lambda2 = lambda;
  }

  function boundsLineStart() {
    boundsStream.point = linePoint;
  }

  function boundsLineEnd() {
    range[0] = lambda0$1, range[1] = lambda1;
    boundsStream.point = boundsPoint;
    p0 = null;
  }

  function boundsRingPoint(lambda, phi) {
    if (p0) {
      var delta = lambda - lambda2;
      deltaSum.add(abs(delta) > 180 ? delta + (delta > 0 ? 360 : -360) : delta);
    } else {
      lambda00$1 = lambda, phi00$1 = phi;
    }
    areaStream.point(lambda, phi);
    linePoint(lambda, phi);
  }

  function boundsRingStart() {
    areaStream.lineStart();
  }

  function boundsRingEnd() {
    boundsRingPoint(lambda00$1, phi00$1);
    areaStream.lineEnd();
    if (abs(deltaSum) > epsilon) lambda0$1 = -(lambda1 = 180);
    range[0] = lambda0$1, range[1] = lambda1;
    p0 = null;
  }

  // Finds the left-right distance between two longitudes.
  // This is almost the same as (lambda1 - lambda0 + 360°) % 360°, except that we want
  // the distance between ±180° to be 360°.
  function angle(lambda0, lambda1) {
    return (lambda1 -= lambda0) < 0 ? lambda1 + 360 : lambda1;
  }

  function rangeCompare(a, b) {
    return a[0] - b[0];
  }

  function rangeContains(range, x) {
    return range[0] <= range[1] ? range[0] <= x && x <= range[1] : x < range[0] || range[1] < x;
  }

  function bounds(feature) {
    var i, n, a, b, merged, deltaMax, delta;

    phi1 = lambda1 = -(lambda0$1 = phi0 = Infinity);
    ranges = [];
    geoStream(feature, boundsStream);

    // First, sort ranges by their minimum longitudes.
    if (n = ranges.length) {
      ranges.sort(rangeCompare);

      // Then, merge any ranges that overlap.
      for (i = 1, a = ranges[0], merged = [a]; i < n; ++i) {
        b = ranges[i];
        if (rangeContains(a, b[0]) || rangeContains(a, b[1])) {
          if (angle(a[0], b[1]) > angle(a[0], a[1])) a[1] = b[1];
          if (angle(b[0], a[1]) > angle(a[0], a[1])) a[0] = b[0];
        } else {
          merged.push(a = b);
        }
      }

      // Finally, find the largest gap between the merged ranges.
      // The final bounding box will be the inverse of this gap.
      for (deltaMax = -Infinity, n = merged.length - 1, i = 0, a = merged[n]; i <= n; a = b, ++i) {
        b = merged[i];
        if ((delta = angle(a[1], b[0])) > deltaMax) deltaMax = delta, lambda0$1 = b[0], lambda1 = a[1];
      }
    }

    ranges = range = null;

    return lambda0$1 === Infinity || phi0 === Infinity
        ? [[NaN, NaN], [NaN, NaN]]
        : [[lambda0$1, phi0], [lambda1, phi1]];
  }

  var W0, W1,
      X0, Y0, Z0,
      X1, Y1, Z1,
      X2, Y2, Z2,
      lambda00$2, phi00$2, // first point
      x0, y0, z0; // previous point

  var centroidStream = {
    sphere: noop,
    point: centroidPoint,
    lineStart: centroidLineStart,
    lineEnd: centroidLineEnd,
    polygonStart: function() {
      centroidStream.lineStart = centroidRingStart;
      centroidStream.lineEnd = centroidRingEnd;
    },
    polygonEnd: function() {
      centroidStream.lineStart = centroidLineStart;
      centroidStream.lineEnd = centroidLineEnd;
    }
  };

  // Arithmetic mean of Cartesian vectors.
  function centroidPoint(lambda, phi) {
    lambda *= radians, phi *= radians;
    var cosPhi = cos(phi);
    centroidPointCartesian(cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi));
  }

  function centroidPointCartesian(x, y, z) {
    ++W0;
    X0 += (x - X0) / W0;
    Y0 += (y - Y0) / W0;
    Z0 += (z - Z0) / W0;
  }

  function centroidLineStart() {
    centroidStream.point = centroidLinePointFirst;
  }

  function centroidLinePointFirst(lambda, phi) {
    lambda *= radians, phi *= radians;
    var cosPhi = cos(phi);
    x0 = cosPhi * cos(lambda);
    y0 = cosPhi * sin(lambda);
    z0 = sin(phi);
    centroidStream.point = centroidLinePoint;
    centroidPointCartesian(x0, y0, z0);
  }

  function centroidLinePoint(lambda, phi) {
    lambda *= radians, phi *= radians;
    var cosPhi = cos(phi),
        x = cosPhi * cos(lambda),
        y = cosPhi * sin(lambda),
        z = sin(phi),
        w = atan2(sqrt((w = y0 * z - z0 * y) * w + (w = z0 * x - x0 * z) * w + (w = x0 * y - y0 * x) * w), x0 * x + y0 * y + z0 * z);
    W1 += w;
    X1 += w * (x0 + (x0 = x));
    Y1 += w * (y0 + (y0 = y));
    Z1 += w * (z0 + (z0 = z));
    centroidPointCartesian(x0, y0, z0);
  }

  function centroidLineEnd() {
    centroidStream.point = centroidPoint;
  }

  // See J. E. Brock, The Inertia Tensor for a Spherical Triangle,
  // J. Applied Mechanics 42, 239 (1975).
  function centroidRingStart() {
    centroidStream.point = centroidRingPointFirst;
  }

  function centroidRingEnd() {
    centroidRingPoint(lambda00$2, phi00$2);
    centroidStream.point = centroidPoint;
  }

  function centroidRingPointFirst(lambda, phi) {
    lambda00$2 = lambda, phi00$2 = phi;
    lambda *= radians, phi *= radians;
    centroidStream.point = centroidRingPoint;
    var cosPhi = cos(phi);
    x0 = cosPhi * cos(lambda);
    y0 = cosPhi * sin(lambda);
    z0 = sin(phi);
    centroidPointCartesian(x0, y0, z0);
  }

  function centroidRingPoint(lambda, phi) {
    lambda *= radians, phi *= radians;
    var cosPhi = cos(phi),
        x = cosPhi * cos(lambda),
        y = cosPhi * sin(lambda),
        z = sin(phi),
        cx = y0 * z - z0 * y,
        cy = z0 * x - x0 * z,
        cz = x0 * y - y0 * x,
        m = sqrt(cx * cx + cy * cy + cz * cz),
        w = asin(m), // line weight = angle
        v = m && -w / m; // area weight multiplier
    X2 += v * cx;
    Y2 += v * cy;
    Z2 += v * cz;
    W1 += w;
    X1 += w * (x0 + (x0 = x));
    Y1 += w * (y0 + (y0 = y));
    Z1 += w * (z0 + (z0 = z));
    centroidPointCartesian(x0, y0, z0);
  }

  function centroid(object) {
    W0 = W1 =
    X0 = Y0 = Z0 =
    X1 = Y1 = Z1 =
    X2 = Y2 = Z2 = 0;
    geoStream(object, centroidStream);

    var x = X2,
        y = Y2,
        z = Z2,
        m = x * x + y * y + z * z;

    // If the area-weighted ccentroid is undefined, fall back to length-weighted ccentroid.
    if (m < epsilon2) {
      x = X1, y = Y1, z = Z1;
      // If the feature has zero length, fall back to arithmetic mean of point vectors.
      if (W1 < epsilon) x = X0, y = Y0, z = Z0;
      m = x * x + y * y + z * z;
      // If the feature still has an undefined ccentroid, then return.
      if (m < epsilon2) return [NaN, NaN];
    }

    return [atan2(y, x) * degrees, asin(z / sqrt(m)) * degrees];
  }

  function constant(x) {
    return function() {
      return x;
    };
  }

  function compose(a, b) {

    function compose(x, y) {
      return x = a(x, y), b(x[0], x[1]);
    }

    if (a.invert && b.invert) compose.invert = function(x, y) {
      return x = b.invert(x, y), x && a.invert(x[0], x[1]);
    };

    return compose;
  }

  function rotationIdentity(lambda, phi) {
    return [lambda > pi ? lambda - tau : lambda < -pi ? lambda + tau : lambda, phi];
  }

  rotationIdentity.invert = rotationIdentity;

  function rotateRadians(deltaLambda, deltaPhi, deltaGamma) {
    return (deltaLambda %= tau) ? (deltaPhi || deltaGamma ? compose(rotationLambda(deltaLambda), rotationPhiGamma(deltaPhi, deltaGamma))
      : rotationLambda(deltaLambda))
      : (deltaPhi || deltaGamma ? rotationPhiGamma(deltaPhi, deltaGamma)
      : rotationIdentity);
  }

  function forwardRotationLambda(deltaLambda) {
    return function(lambda, phi) {
      return lambda += deltaLambda, [lambda > pi ? lambda - tau : lambda < -pi ? lambda + tau : lambda, phi];
    };
  }

  function rotationLambda(deltaLambda) {
    var rotation = forwardRotationLambda(deltaLambda);
    rotation.invert = forwardRotationLambda(-deltaLambda);
    return rotation;
  }

  function rotationPhiGamma(deltaPhi, deltaGamma) {
    var cosDeltaPhi = cos(deltaPhi),
        sinDeltaPhi = sin(deltaPhi),
        cosDeltaGamma = cos(deltaGamma),
        sinDeltaGamma = sin(deltaGamma);

    function rotation(lambda, phi) {
      var cosPhi = cos(phi),
          x = cos(lambda) * cosPhi,
          y = sin(lambda) * cosPhi,
          z = sin(phi),
          k = z * cosDeltaPhi + x * sinDeltaPhi;
      return [
        atan2(y * cosDeltaGamma - k * sinDeltaGamma, x * cosDeltaPhi - z * sinDeltaPhi),
        asin(k * cosDeltaGamma + y * sinDeltaGamma)
      ];
    }

    rotation.invert = function(lambda, phi) {
      var cosPhi = cos(phi),
          x = cos(lambda) * cosPhi,
          y = sin(lambda) * cosPhi,
          z = sin(phi),
          k = z * cosDeltaGamma - y * sinDeltaGamma;
      return [
        atan2(y * cosDeltaGamma + z * sinDeltaGamma, x * cosDeltaPhi + k * sinDeltaPhi),
        asin(k * cosDeltaPhi - x * sinDeltaPhi)
      ];
    };

    return rotation;
  }

  function rotation(rotate) {
    rotate = rotateRadians(rotate[0] * radians, rotate[1] * radians, rotate.length > 2 ? rotate[2] * radians : 0);

    function forward(coordinates) {
      coordinates = rotate(coordinates[0] * radians, coordinates[1] * radians);
      return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
    }

    forward.invert = function(coordinates) {
      coordinates = rotate.invert(coordinates[0] * radians, coordinates[1] * radians);
      return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
    };

    return forward;
  }

  // Generates a circle centered at [0°, 0°], with a given radius and precision.
  function circleStream(stream, radius, delta, direction, t0, t1) {
    if (!delta) return;
    var cosRadius = cos(radius),
        sinRadius = sin(radius),
        step = direction * delta;
    if (t0 == null) {
      t0 = radius + direction * tau;
      t1 = radius - step / 2;
    } else {
      t0 = circleRadius(cosRadius, t0);
      t1 = circleRadius(cosRadius, t1);
      if (direction > 0 ? t0 < t1 : t0 > t1) t0 += direction * tau;
    }
    for (var point, t = t0; direction > 0 ? t > t1 : t < t1; t -= step) {
      point = spherical([cosRadius, -sinRadius * cos(t), -sinRadius * sin(t)]);
      stream.point(point[0], point[1]);
    }
  }

  // Returns the signed angle of a cartesian point relative to [cosRadius, 0, 0].
  function circleRadius(cosRadius, point) {
    point = cartesian(point), point[0] -= cosRadius;
    cartesianNormalizeInPlace(point);
    var radius = acos(-point[1]);
    return ((-point[2] < 0 ? -radius : radius) + tau - epsilon) % tau;
  }

  function geoCircle() {
    var center = constant([0, 0]),
        radius = constant(90),
        precision = constant(6),
        ring,
        rotate,
        stream = {point: point};

    function point(x, y) {
      ring.push(x = rotate(x, y));
      x[0] *= degrees, x[1] *= degrees;
    }

    function circle() {
      var c = center.apply(this, arguments),
          r = radius.apply(this, arguments) * radians,
          p = precision.apply(this, arguments) * radians;
      ring = [];
      rotate = rotateRadians(-c[0] * radians, -c[1] * radians, 0).invert;
      circleStream(stream, r, p, 1);
      c = {type: "Polygon", coordinates: [ring]};
      ring = rotate = null;
      return c;
    }

    circle.center = function(_) {
      return arguments.length ? (center = typeof _ === "function" ? _ : constant([+_[0], +_[1]]), circle) : center;
    };

    circle.radius = function(_) {
      return arguments.length ? (radius = typeof _ === "function" ? _ : constant(+_), circle) : radius;
    };

    circle.precision = function(_) {
      return arguments.length ? (precision = typeof _ === "function" ? _ : constant(+_), circle) : precision;
    };

    return circle;
  }

  function clipBuffer() {
    var lines = [],
        line;
    return {
      point: function(x, y) {
        line.push([x, y]);
      },
      lineStart: function() {
        lines.push(line = []);
      },
      lineEnd: noop,
      rejoin: function() {
        if (lines.length > 1) lines.push(lines.pop().concat(lines.shift()));
      },
      result: function() {
        var result = lines;
        lines = [];
        line = null;
        return result;
      }
    };
  }

  function pointEqual(a, b) {
    return abs(a[0] - b[0]) < epsilon && abs(a[1] - b[1]) < epsilon;
  }

  function Intersection(point, points, other, entry) {
    this.x = point;
    this.z = points;
    this.o = other; // another intersection
    this.e = entry; // is an entry?
    this.v = false; // visited
    this.n = this.p = null; // next & previous
  }

  // A generalized polygon clipping algorithm: given a polygon that has been cut
  // into its visible line segments, and rejoins the segments by interpolating
  // along the clip edge.
  function clipRejoin(segments, compareIntersection, startInside, interpolate, stream) {
    var subject = [],
        clip = [],
        i,
        n;

    segments.forEach(function(segment) {
      if ((n = segment.length - 1) <= 0) return;
      var n, p0 = segment[0], p1 = segment[n], x;

      // If the first and last points of a segment are coincident, then treat as a
      // closed ring. TODO if all rings are closed, then the winding order of the
      // exterior ring should be checked.
      if (pointEqual(p0, p1)) {
        stream.lineStart();
        for (i = 0; i < n; ++i) stream.point((p0 = segment[i])[0], p0[1]);
        stream.lineEnd();
        return;
      }

      subject.push(x = new Intersection(p0, segment, null, true));
      clip.push(x.o = new Intersection(p0, null, x, false));
      subject.push(x = new Intersection(p1, segment, null, false));
      clip.push(x.o = new Intersection(p1, null, x, true));
    });

    if (!subject.length) return;

    clip.sort(compareIntersection);
    link(subject);
    link(clip);

    for (i = 0, n = clip.length; i < n; ++i) {
      clip[i].e = startInside = !startInside;
    }

    var start = subject[0],
        points,
        point;

    while (1) {
      // Find first unvisited intersection.
      var current = start,
          isSubject = true;
      while (current.v) if ((current = current.n) === start) return;
      points = current.z;
      stream.lineStart();
      do {
        current.v = current.o.v = true;
        if (current.e) {
          if (isSubject) {
            for (i = 0, n = points.length; i < n; ++i) stream.point((point = points[i])[0], point[1]);
          } else {
            interpolate(current.x, current.n.x, 1, stream);
          }
          current = current.n;
        } else {
          if (isSubject) {
            points = current.p.z;
            for (i = points.length - 1; i >= 0; --i) stream.point((point = points[i])[0], point[1]);
          } else {
            interpolate(current.x, current.p.x, -1, stream);
          }
          current = current.p;
        }
        current = current.o;
        points = current.z;
        isSubject = !isSubject;
      } while (!current.v);
      stream.lineEnd();
    }
  }

  function link(array) {
    if (!(n = array.length)) return;
    var n,
        i = 0,
        a = array[0],
        b;
    while (++i < n) {
      a.n = b = array[i];
      b.p = a;
      a = b;
    }
    a.n = b = array[0];
    b.p = a;
  }

  var sum = adder();

  function polygonContains(polygon, point) {
    var lambda = point[0],
        phi = point[1],
        sinPhi = sin(phi),
        normal = [sin(lambda), -cos(lambda), 0],
        angle = 0,
        winding = 0;

    sum.reset();

    if (sinPhi === 1) phi = halfPi + epsilon;
    else if (sinPhi === -1) phi = -halfPi - epsilon;

    for (var i = 0, n = polygon.length; i < n; ++i) {
      if (!(m = (ring = polygon[i]).length)) continue;
      var ring,
          m,
          point0 = ring[m - 1],
          lambda0 = point0[0],
          phi0 = point0[1] / 2 + quarterPi,
          sinPhi0 = sin(phi0),
          cosPhi0 = cos(phi0);

      for (var j = 0; j < m; ++j, lambda0 = lambda1, sinPhi0 = sinPhi1, cosPhi0 = cosPhi1, point0 = point1) {
        var point1 = ring[j],
            lambda1 = point1[0],
            phi1 = point1[1] / 2 + quarterPi,
            sinPhi1 = sin(phi1),
            cosPhi1 = cos(phi1),
            delta = lambda1 - lambda0,
            sign$$1 = delta >= 0 ? 1 : -1,
            absDelta = sign$$1 * delta,
            antimeridian = absDelta > pi,
            k = sinPhi0 * sinPhi1;

        sum.add(atan2(k * sign$$1 * sin(absDelta), cosPhi0 * cosPhi1 + k * cos(absDelta)));
        angle += antimeridian ? delta + sign$$1 * tau : delta;

        // Are the longitudes either side of the point’s meridian (lambda),
        // and are the latitudes smaller than the parallel (phi)?
        if (antimeridian ^ lambda0 >= lambda ^ lambda1 >= lambda) {
          var arc = cartesianCross(cartesian(point0), cartesian(point1));
          cartesianNormalizeInPlace(arc);
          var intersection = cartesianCross(normal, arc);
          cartesianNormalizeInPlace(intersection);
          var phiArc = (antimeridian ^ delta >= 0 ? -1 : 1) * asin(intersection[2]);
          if (phi > phiArc || phi === phiArc && (arc[0] || arc[1])) {
            winding += antimeridian ^ delta >= 0 ? 1 : -1;
          }
        }
      }
    }

    // First, determine whether the South pole is inside or outside:
    //
    // It is inside if:
    // * the polygon winds around it in a clockwise direction.
    // * the polygon does not (cumulatively) wind around it, but has a negative
    //   (counter-clockwise) area.
    //
    // Second, count the (signed) number of times a segment crosses a lambda
    // from the point to the South pole.  If it is zero, then the point is the
    // same side as the South pole.

    return (angle < -epsilon || angle < epsilon && sum < -epsilon) ^ (winding & 1);
  }

  function ascending(a, b) {
    return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
  }

  function bisector(compare) {
    if (compare.length === 1) compare = ascendingComparator(compare);
    return {
      left: function(a, x, lo, hi) {
        if (lo == null) lo = 0;
        if (hi == null) hi = a.length;
        while (lo < hi) {
          var mid = lo + hi >>> 1;
          if (compare(a[mid], x) < 0) lo = mid + 1;
          else hi = mid;
        }
        return lo;
      },
      right: function(a, x, lo, hi) {
        if (lo == null) lo = 0;
        if (hi == null) hi = a.length;
        while (lo < hi) {
          var mid = lo + hi >>> 1;
          if (compare(a[mid], x) > 0) hi = mid;
          else lo = mid + 1;
        }
        return lo;
      }
    };
  }

  function ascendingComparator(f) {
    return function(d, x) {
      return ascending(f(d), x);
    };
  }

  var ascendingBisect = bisector(ascending);

  function range$1(start, stop, step) {
    start = +start, stop = +stop, step = (n = arguments.length) < 2 ? (stop = start, start = 0, 1) : n < 3 ? 1 : +step;

    var i = -1,
        n = Math.max(0, Math.ceil((stop - start) / step)) | 0,
        range = new Array(n);

    while (++i < n) {
      range[i] = start + i * step;
    }

    return range;
  }

  function merge(arrays) {
    var n = arrays.length,
        m,
        i = -1,
        j = 0,
        merged,
        array;

    while (++i < n) j += arrays[i].length;
    merged = new Array(j);

    while (--n >= 0) {
      array = arrays[n];
      m = array.length;
      while (--m >= 0) {
        merged[--j] = array[m];
      }
    }

    return merged;
  }

  function clip(pointVisible, clipLine, interpolate, start) {
    return function(sink) {
      var line = clipLine(sink),
          ringBuffer = clipBuffer(),
          ringSink = clipLine(ringBuffer),
          polygonStarted = false,
          polygon,
          segments,
          ring;

      var clip = {
        point: point,
        lineStart: lineStart,
        lineEnd: lineEnd,
        polygonStart: function() {
          clip.point = pointRing;
          clip.lineStart = ringStart;
          clip.lineEnd = ringEnd;
          segments = [];
          polygon = [];
        },
        polygonEnd: function() {
          clip.point = point;
          clip.lineStart = lineStart;
          clip.lineEnd = lineEnd;
          segments = merge(segments);
          var startInside = polygonContains(polygon, start);
          if (segments.length) {
            if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
            clipRejoin(segments, compareIntersection, startInside, interpolate, sink);
          } else if (startInside) {
            if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
            sink.lineStart();
            interpolate(null, null, 1, sink);
            sink.lineEnd();
          }
          if (polygonStarted) sink.polygonEnd(), polygonStarted = false;
          segments = polygon = null;
        },
        sphere: function() {
          sink.polygonStart();
          sink.lineStart();
          interpolate(null, null, 1, sink);
          sink.lineEnd();
          sink.polygonEnd();
        }
      };

      function point(lambda, phi) {
        if (pointVisible(lambda, phi)) sink.point(lambda, phi);
      }

      function pointLine(lambda, phi) {
        line.point(lambda, phi);
      }

      function lineStart() {
        clip.point = pointLine;
        line.lineStart();
      }

      function lineEnd() {
        clip.point = point;
        line.lineEnd();
      }

      function pointRing(lambda, phi) {
        ring.push([lambda, phi]);
        ringSink.point(lambda, phi);
      }

      function ringStart() {
        ringSink.lineStart();
        ring = [];
      }

      function ringEnd() {
        pointRing(ring[0][0], ring[0][1]);
        ringSink.lineEnd();

        var clean = ringSink.clean(),
            ringSegments = ringBuffer.result(),
            i, n = ringSegments.length, m,
            segment,
            point;

        ring.pop();
        polygon.push(ring);
        ring = null;

        if (!n) return;

        // No intersections.
        if (clean & 1) {
          segment = ringSegments[0];
          if ((m = segment.length - 1) > 0) {
            if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
            sink.lineStart();
            for (i = 0; i < m; ++i) sink.point((point = segment[i])[0], point[1]);
            sink.lineEnd();
          }
          return;
        }

        // Rejoin connected segments.
        // TODO reuse ringBuffer.rejoin()?
        if (n > 1 && clean & 2) ringSegments.push(ringSegments.pop().concat(ringSegments.shift()));

        segments.push(ringSegments.filter(validSegment));
      }

      return clip;
    };
  }

  function validSegment(segment) {
    return segment.length > 1;
  }

  // Intersections are sorted along the clip edge. For both antimeridian cutting
  // and circle clipping, the same comparison is used.
  function compareIntersection(a, b) {
    return ((a = a.x)[0] < 0 ? a[1] - halfPi - epsilon : halfPi - a[1])
         - ((b = b.x)[0] < 0 ? b[1] - halfPi - epsilon : halfPi - b[1]);
  }

  var clipAntimeridian = clip(
    function() { return true; },
    clipAntimeridianLine,
    clipAntimeridianInterpolate,
    [-pi, -halfPi]
  );

  // Takes a line and cuts into visible segments. Return values: 0 - there were
  // intersections or the line was empty; 1 - no intersections; 2 - there were
  // intersections, and the first and last segments should be rejoined.
  function clipAntimeridianLine(stream) {
    var lambda0 = NaN,
        phi0 = NaN,
        sign0 = NaN,
        clean; // no intersections

    return {
      lineStart: function() {
        stream.lineStart();
        clean = 1;
      },
      point: function(lambda1, phi1) {
        var sign1 = lambda1 > 0 ? pi : -pi,
            delta = abs(lambda1 - lambda0);
        if (abs(delta - pi) < epsilon) { // line crosses a pole
          stream.point(lambda0, phi0 = (phi0 + phi1) / 2 > 0 ? halfPi : -halfPi);
          stream.point(sign0, phi0);
          stream.lineEnd();
          stream.lineStart();
          stream.point(sign1, phi0);
          stream.point(lambda1, phi0);
          clean = 0;
        } else if (sign0 !== sign1 && delta >= pi) { // line crosses antimeridian
          if (abs(lambda0 - sign0) < epsilon) lambda0 -= sign0 * epsilon; // handle degeneracies
          if (abs(lambda1 - sign1) < epsilon) lambda1 -= sign1 * epsilon;
          phi0 = clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1);
          stream.point(sign0, phi0);
          stream.lineEnd();
          stream.lineStart();
          stream.point(sign1, phi0);
          clean = 0;
        }
        stream.point(lambda0 = lambda1, phi0 = phi1);
        sign0 = sign1;
      },
      lineEnd: function() {
        stream.lineEnd();
        lambda0 = phi0 = NaN;
      },
      clean: function() {
        return 2 - clean; // if intersections, rejoin first and last segments
      }
    };
  }

  function clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1) {
    var cosPhi0,
        cosPhi1,
        sinLambda0Lambda1 = sin(lambda0 - lambda1);
    return abs(sinLambda0Lambda1) > epsilon
        ? atan((sin(phi0) * (cosPhi1 = cos(phi1)) * sin(lambda1)
            - sin(phi1) * (cosPhi0 = cos(phi0)) * sin(lambda0))
            / (cosPhi0 * cosPhi1 * sinLambda0Lambda1))
        : (phi0 + phi1) / 2;
  }

  function clipAntimeridianInterpolate(from, to, direction, stream) {
    var phi;
    if (from == null) {
      phi = direction * halfPi;
      stream.point(-pi, phi);
      stream.point(0, phi);
      stream.point(pi, phi);
      stream.point(pi, 0);
      stream.point(pi, -phi);
      stream.point(0, -phi);
      stream.point(-pi, -phi);
      stream.point(-pi, 0);
      stream.point(-pi, phi);
    } else if (abs(from[0] - to[0]) > epsilon) {
      var lambda = from[0] < to[0] ? pi : -pi;
      phi = direction * lambda / 2;
      stream.point(-lambda, phi);
      stream.point(0, phi);
      stream.point(lambda, phi);
    } else {
      stream.point(to[0], to[1]);
    }
  }

  function clipCircle(radius) {
    var cr = cos(radius),
        delta = 6 * radians,
        smallRadius = cr > 0,
        notHemisphere = abs(cr) > epsilon; // TODO optimise for this common case

    function interpolate(from, to, direction, stream) {
      circleStream(stream, radius, delta, direction, from, to);
    }

    function visible(lambda, phi) {
      return cos(lambda) * cos(phi) > cr;
    }

    // Takes a line and cuts into visible segments. Return values used for polygon
    // clipping: 0 - there were intersections or the line was empty; 1 - no
    // intersections 2 - there were intersections, and the first and last segments
    // should be rejoined.
    function clipLine(stream) {
      var point0, // previous point
          c0, // code for previous point
          v0, // visibility of previous point
          v00, // visibility of first point
          clean; // no intersections
      return {
        lineStart: function() {
          v00 = v0 = false;
          clean = 1;
        },
        point: function(lambda, phi) {
          var point1 = [lambda, phi],
              point2,
              v = visible(lambda, phi),
              c = smallRadius
                ? v ? 0 : code(lambda, phi)
                : v ? code(lambda + (lambda < 0 ? pi : -pi), phi) : 0;
          if (!point0 && (v00 = v0 = v)) stream.lineStart();
          // Handle degeneracies.
          // TODO ignore if not clipping polygons.
          if (v !== v0) {
            point2 = intersect(point0, point1);
            if (!point2 || pointEqual(point0, point2) || pointEqual(point1, point2)) {
              point1[0] += epsilon;
              point1[1] += epsilon;
              v = visible(point1[0], point1[1]);
            }
          }
          if (v !== v0) {
            clean = 0;
            if (v) {
              // outside going in
              stream.lineStart();
              point2 = intersect(point1, point0);
              stream.point(point2[0], point2[1]);
            } else {
              // inside going out
              point2 = intersect(point0, point1);
              stream.point(point2[0], point2[1]);
              stream.lineEnd();
            }
            point0 = point2;
          } else if (notHemisphere && point0 && smallRadius ^ v) {
            var t;
            // If the codes for two points are different, or are both zero,
            // and there this segment intersects with the small circle.
            if (!(c & c0) && (t = intersect(point1, point0, true))) {
              clean = 0;
              if (smallRadius) {
                stream.lineStart();
                stream.point(t[0][0], t[0][1]);
                stream.point(t[1][0], t[1][1]);
                stream.lineEnd();
              } else {
                stream.point(t[1][0], t[1][1]);
                stream.lineEnd();
                stream.lineStart();
                stream.point(t[0][0], t[0][1]);
              }
            }
          }
          if (v && (!point0 || !pointEqual(point0, point1))) {
            stream.point(point1[0], point1[1]);
          }
          point0 = point1, v0 = v, c0 = c;
        },
        lineEnd: function() {
          if (v0) stream.lineEnd();
          point0 = null;
        },
        // Rejoin first and last segments if there were intersections and the first
        // and last points were visible.
        clean: function() {
          return clean | ((v00 && v0) << 1);
        }
      };
    }

    // Intersects the great circle between a and b with the clip circle.
    function intersect(a, b, two) {
      var pa = cartesian(a),
          pb = cartesian(b);

      // We have two planes, n1.p = d1 and n2.p = d2.
      // Find intersection line p(t) = c1 n1 + c2 n2 + t (n1 ⨯ n2).
      var n1 = [1, 0, 0], // normal
          n2 = cartesianCross(pa, pb),
          n2n2 = cartesianDot(n2, n2),
          n1n2 = n2[0], // cartesianDot(n1, n2),
          determinant = n2n2 - n1n2 * n1n2;

      // Two polar points.
      if (!determinant) return !two && a;

      var c1 =  cr * n2n2 / determinant,
          c2 = -cr * n1n2 / determinant,
          n1xn2 = cartesianCross(n1, n2),
          A = cartesianScale(n1, c1),
          B = cartesianScale(n2, c2);
      cartesianAddInPlace(A, B);

      // Solve |p(t)|^2 = 1.
      var u = n1xn2,
          w = cartesianDot(A, u),
          uu = cartesianDot(u, u),
          t2 = w * w - uu * (cartesianDot(A, A) - 1);

      if (t2 < 0) return;

      var t = sqrt(t2),
          q = cartesianScale(u, (-w - t) / uu);
      cartesianAddInPlace(q, A);
      q = spherical(q);

      if (!two) return q;

      // Two intersection points.
      var lambda0 = a[0],
          lambda1 = b[0],
          phi0 = a[1],
          phi1 = b[1],
          z;

      if (lambda1 < lambda0) z = lambda0, lambda0 = lambda1, lambda1 = z;

      var delta = lambda1 - lambda0,
          polar = abs(delta - pi) < epsilon,
          meridian = polar || delta < epsilon;

      if (!polar && phi1 < phi0) z = phi0, phi0 = phi1, phi1 = z;

      // Check that the first point is between a and b.
      if (meridian
          ? polar
            ? phi0 + phi1 > 0 ^ q[1] < (abs(q[0] - lambda0) < epsilon ? phi0 : phi1)
            : phi0 <= q[1] && q[1] <= phi1
          : delta > pi ^ (lambda0 <= q[0] && q[0] <= lambda1)) {
        var q1 = cartesianScale(u, (-w + t) / uu);
        cartesianAddInPlace(q1, A);
        return [q, spherical(q1)];
      }
    }

    // Generates a 4-bit vector representing the location of a point relative to
    // the small circle's bounding box.
    function code(lambda, phi) {
      var r = smallRadius ? radius : pi - radius,
          code = 0;
      if (lambda < -r) code |= 1; // left
      else if (lambda > r) code |= 2; // right
      if (phi < -r) code |= 4; // below
      else if (phi > r) code |= 8; // above
      return code;
    }

    return clip(visible, clipLine, interpolate, smallRadius ? [0, -radius] : [-pi, radius - pi]);
  }

  function clipLine(a, b, x0, y0, x1, y1) {
    var ax = a[0],
        ay = a[1],
        bx = b[0],
        by = b[1],
        t0 = 0,
        t1 = 1,
        dx = bx - ax,
        dy = by - ay,
        r;

    r = x0 - ax;
    if (!dx && r > 0) return;
    r /= dx;
    if (dx < 0) {
      if (r < t0) return;
      if (r < t1) t1 = r;
    } else if (dx > 0) {
      if (r > t1) return;
      if (r > t0) t0 = r;
    }

    r = x1 - ax;
    if (!dx && r < 0) return;
    r /= dx;
    if (dx < 0) {
      if (r > t1) return;
      if (r > t0) t0 = r;
    } else if (dx > 0) {
      if (r < t0) return;
      if (r < t1) t1 = r;
    }

    r = y0 - ay;
    if (!dy && r > 0) return;
    r /= dy;
    if (dy < 0) {
      if (r < t0) return;
      if (r < t1) t1 = r;
    } else if (dy > 0) {
      if (r > t1) return;
      if (r > t0) t0 = r;
    }

    r = y1 - ay;
    if (!dy && r < 0) return;
    r /= dy;
    if (dy < 0) {
      if (r > t1) return;
      if (r > t0) t0 = r;
    } else if (dy > 0) {
      if (r < t0) return;
      if (r < t1) t1 = r;
    }

    if (t0 > 0) a[0] = ax + t0 * dx, a[1] = ay + t0 * dy;
    if (t1 < 1) b[0] = ax + t1 * dx, b[1] = ay + t1 * dy;
    return true;
  }

  var clipMax = 1e9, clipMin = -clipMax;

  // TODO Use d3-polygon’s polygonContains here for the ring check?
  // TODO Eliminate duplicate buffering in clipBuffer and polygon.push?

  function clipRectangle(x0, y0, x1, y1) {

    function visible(x, y) {
      return x0 <= x && x <= x1 && y0 <= y && y <= y1;
    }

    function interpolate(from, to, direction, stream) {
      var a = 0, a1 = 0;
      if (from == null
          || (a = corner(from, direction)) !== (a1 = corner(to, direction))
          || comparePoint(from, to) < 0 ^ direction > 0) {
        do stream.point(a === 0 || a === 3 ? x0 : x1, a > 1 ? y1 : y0);
        while ((a = (a + direction + 4) % 4) !== a1);
      } else {
        stream.point(to[0], to[1]);
      }
    }

    function corner(p, direction) {
      return abs(p[0] - x0) < epsilon ? direction > 0 ? 0 : 3
          : abs(p[0] - x1) < epsilon ? direction > 0 ? 2 : 1
          : abs(p[1] - y0) < epsilon ? direction > 0 ? 1 : 0
          : direction > 0 ? 3 : 2; // abs(p[1] - y1) < epsilon
    }

    function compareIntersection(a, b) {
      return comparePoint(a.x, b.x);
    }

    function comparePoint(a, b) {
      var ca = corner(a, 1),
          cb = corner(b, 1);
      return ca !== cb ? ca - cb
          : ca === 0 ? b[1] - a[1]
          : ca === 1 ? a[0] - b[0]
          : ca === 2 ? a[1] - b[1]
          : b[0] - a[0];
    }

    return function(stream) {
      var activeStream = stream,
          bufferStream = clipBuffer(),
          segments,
          polygon,
          ring,
          x__, y__, v__, // first point
          x_, y_, v_, // previous point
          first,
          clean;

      var clipStream = {
        point: point,
        lineStart: lineStart,
        lineEnd: lineEnd,
        polygonStart: polygonStart,
        polygonEnd: polygonEnd
      };

      function point(x, y) {
        if (visible(x, y)) activeStream.point(x, y);
      }

      function polygonInside() {
        var winding = 0;

        for (var i = 0, n = polygon.length; i < n; ++i) {
          for (var ring = polygon[i], j = 1, m = ring.length, point = ring[0], a0, a1, b0 = point[0], b1 = point[1]; j < m; ++j) {
            a0 = b0, a1 = b1, point = ring[j], b0 = point[0], b1 = point[1];
            if (a1 <= y1) { if (b1 > y1 && (b0 - a0) * (y1 - a1) > (b1 - a1) * (x0 - a0)) ++winding; }
            else { if (b1 <= y1 && (b0 - a0) * (y1 - a1) < (b1 - a1) * (x0 - a0)) --winding; }
          }
        }

        return winding;
      }

      // Buffer geometry within a polygon and then clip it en masse.
      function polygonStart() {
        activeStream = bufferStream, segments = [], polygon = [], clean = true;
      }

      function polygonEnd() {
        var startInside = polygonInside(),
            cleanInside = clean && startInside,
            visible = (segments = merge(segments)).length;
        if (cleanInside || visible) {
          stream.polygonStart();
          if (cleanInside) {
            stream.lineStart();
            interpolate(null, null, 1, stream);
            stream.lineEnd();
          }
          if (visible) {
            clipRejoin(segments, compareIntersection, startInside, interpolate, stream);
          }
          stream.polygonEnd();
        }
        activeStream = stream, segments = polygon = ring = null;
      }

      function lineStart() {
        clipStream.point = linePoint;
        if (polygon) polygon.push(ring = []);
        first = true;
        v_ = false;
        x_ = y_ = NaN;
      }

      // TODO rather than special-case polygons, simply handle them separately.
      // Ideally, coincident intersection points should be jittered to avoid
      // clipping issues.
      function lineEnd() {
        if (segments) {
          linePoint(x__, y__);
          if (v__ && v_) bufferStream.rejoin();
          segments.push(bufferStream.result());
        }
        clipStream.point = point;
        if (v_) activeStream.lineEnd();
      }

      function linePoint(x, y) {
        var v = visible(x, y);
        if (polygon) ring.push([x, y]);
        if (first) {
          x__ = x, y__ = y, v__ = v;
          first = false;
          if (v) {
            activeStream.lineStart();
            activeStream.point(x, y);
          }
        } else {
          if (v && v_) activeStream.point(x, y);
          else {
            var a = [x_ = Math.max(clipMin, Math.min(clipMax, x_)), y_ = Math.max(clipMin, Math.min(clipMax, y_))],
                b = [x = Math.max(clipMin, Math.min(clipMax, x)), y = Math.max(clipMin, Math.min(clipMax, y))];
            if (clipLine(a, b, x0, y0, x1, y1)) {
              if (!v_) {
                activeStream.lineStart();
                activeStream.point(a[0], a[1]);
              }
              activeStream.point(b[0], b[1]);
              if (!v) activeStream.lineEnd();
              clean = false;
            } else if (v) {
              activeStream.lineStart();
              activeStream.point(x, y);
              clean = false;
            }
          }
        }
        x_ = x, y_ = y, v_ = v;
      }

      return clipStream;
    };
  }

  var lengthSum = adder();

  function interpolate(a, b) {
    var x0 = a[0] * radians,
        y0 = a[1] * radians,
        x1 = b[0] * radians,
        y1 = b[1] * radians,
        cy0 = cos(y0),
        sy0 = sin(y0),
        cy1 = cos(y1),
        sy1 = sin(y1),
        kx0 = cy0 * cos(x0),
        ky0 = cy0 * sin(x0),
        kx1 = cy1 * cos(x1),
        ky1 = cy1 * sin(x1),
        d = 2 * asin(sqrt(haversin(y1 - y0) + cy0 * cy1 * haversin(x1 - x0))),
        k = sin(d);

    var interpolate = d ? function(t) {
      var B = sin(t *= d) / k,
          A = sin(d - t) / k,
          x = A * kx0 + B * kx1,
          y = A * ky0 + B * ky1,
          z = A * sy0 + B * sy1;
      return [
        atan2(y, x) * degrees,
        atan2(z, sqrt(x * x + y * y)) * degrees
      ];
    } : function() {
      return [x0 * degrees, y0 * degrees];
    };

    interpolate.distance = d;

    return interpolate;
  }

  function identity$1(x) {
    return x;
  }

  var areaSum$1 = adder(),
      areaRingSum$1 = adder();

  var x0$2 = Infinity,
      y0$2 = x0$2,
      x1 = -x0$2,
      y1 = x1;

  var boundsStream$1 = {
    point: boundsPoint$1,
    lineStart: noop,
    lineEnd: noop,
    polygonStart: noop,
    polygonEnd: noop,
    result: function() {
      var bounds = [[x0$2, y0$2], [x1, y1]];
      x1 = y1 = -(y0$2 = x0$2 = Infinity);
      return bounds;
    }
  };

  function boundsPoint$1(x, y) {
    if (x < x0$2) x0$2 = x;
    if (x > x1) x1 = x;
    if (y < y0$2) y0$2 = y;
    if (y > y1) y1 = y;
  }

  var lengthSum$1 = adder();

  function transformer(methods) {
    return function(stream) {
      var s = new TransformStream;
      for (var key in methods) s[key] = methods[key];
      s.stream = stream;
      return s;
    };
  }

  function TransformStream() {}

  TransformStream.prototype = {
    constructor: TransformStream,
    point: function(x, y) { this.stream.point(x, y); },
    sphere: function() { this.stream.sphere(); },
    lineStart: function() { this.stream.lineStart(); },
    lineEnd: function() { this.stream.lineEnd(); },
    polygonStart: function() { this.stream.polygonStart(); },
    polygonEnd: function() { this.stream.polygonEnd(); }
  };

  function fit(projection, fitBounds, object) {
    var clip = projection.clipExtent && projection.clipExtent();
    projection.scale(150).translate([0, 0]);
    if (clip != null) projection.clipExtent(null);
    geoStream(object, projection.stream(boundsStream$1));
    fitBounds(boundsStream$1.result());
    if (clip != null) projection.clipExtent(clip);
    return projection;
  }

  function fitExtent(projection, extent, object) {
    return fit(projection, function(b) {
      var w = extent[1][0] - extent[0][0],
          h = extent[1][1] - extent[0][1],
          k = Math.min(w / (b[1][0] - b[0][0]), h / (b[1][1] - b[0][1])),
          x = +extent[0][0] + (w - k * (b[1][0] + b[0][0])) / 2,
          y = +extent[0][1] + (h - k * (b[1][1] + b[0][1])) / 2;
      projection.scale(150 * k).translate([x, y]);
    }, object);
  }

  function fitSize(projection, size, object) {
    return fitExtent(projection, [[0, 0], size], object);
  }

  function fitWidth(projection, width, object) {
    return fit(projection, function(b) {
      var w = +width,
          k = w / (b[1][0] - b[0][0]),
          x = (w - k * (b[1][0] + b[0][0])) / 2,
          y = -k * b[0][1];
      projection.scale(150 * k).translate([x, y]);
    }, object);
  }

  function fitHeight(projection, height, object) {
    return fit(projection, function(b) {
      var h = +height,
          k = h / (b[1][1] - b[0][1]),
          x = -k * b[0][0],
          y = (h - k * (b[1][1] + b[0][1])) / 2;
      projection.scale(150 * k).translate([x, y]);
    }, object);
  }

  var maxDepth = 16, // maximum depth of subdivision
      cosMinDistance = cos(30 * radians); // cos(minimum angular distance)

  function resample(project, delta2) {
    return +delta2 ? resample$1(project, delta2) : resampleNone(project);
  }

  function resampleNone(project) {
    return transformer({
      point: function(x, y) {
        x = project(x, y);
        this.stream.point(x[0], x[1]);
      }
    });
  }

  function resample$1(project, delta2) {

    function resampleLineTo(x0, y0, lambda0, a0, b0, c0, x1, y1, lambda1, a1, b1, c1, depth, stream) {
      var dx = x1 - x0,
          dy = y1 - y0,
          d2 = dx * dx + dy * dy;
      if (d2 > 4 * delta2 && depth--) {
        var a = a0 + a1,
            b = b0 + b1,
            c = c0 + c1,
            m = sqrt(a * a + b * b + c * c),
            phi2 = asin(c /= m),
            lambda2 = abs(abs(c) - 1) < epsilon || abs(lambda0 - lambda1) < epsilon ? (lambda0 + lambda1) / 2 : atan2(b, a),
            p = project(lambda2, phi2),
            x2 = p[0],
            y2 = p[1],
            dx2 = x2 - x0,
            dy2 = y2 - y0,
            dz = dy * dx2 - dx * dy2;
        if (dz * dz / d2 > delta2 // perpendicular projected distance
            || abs((dx * dx2 + dy * dy2) / d2 - 0.5) > 0.3 // midpoint close to an end
            || a0 * a1 + b0 * b1 + c0 * c1 < cosMinDistance) { // angular distance
          resampleLineTo(x0, y0, lambda0, a0, b0, c0, x2, y2, lambda2, a /= m, b /= m, c, depth, stream);
          stream.point(x2, y2);
          resampleLineTo(x2, y2, lambda2, a, b, c, x1, y1, lambda1, a1, b1, c1, depth, stream);
        }
      }
    }
    return function(stream) {
      var lambda00, x00, y00, a00, b00, c00, // first point
          lambda0, x0, y0, a0, b0, c0; // previous point

      var resampleStream = {
        point: point,
        lineStart: lineStart,
        lineEnd: lineEnd,
        polygonStart: function() { stream.polygonStart(); resampleStream.lineStart = ringStart; },
        polygonEnd: function() { stream.polygonEnd(); resampleStream.lineStart = lineStart; }
      };

      function point(x, y) {
        x = project(x, y);
        stream.point(x[0], x[1]);
      }

      function lineStart() {
        x0 = NaN;
        resampleStream.point = linePoint;
        stream.lineStart();
      }

      function linePoint(lambda, phi) {
        var c = cartesian([lambda, phi]), p = project(lambda, phi);
        resampleLineTo(x0, y0, lambda0, a0, b0, c0, x0 = p[0], y0 = p[1], lambda0 = lambda, a0 = c[0], b0 = c[1], c0 = c[2], maxDepth, stream);
        stream.point(x0, y0);
      }

      function lineEnd() {
        resampleStream.point = point;
        stream.lineEnd();
      }

      function ringStart() {
        lineStart();
        resampleStream.point = ringPoint;
        resampleStream.lineEnd = ringEnd;
      }

      function ringPoint(lambda, phi) {
        linePoint(lambda00 = lambda, phi), x00 = x0, y00 = y0, a00 = a0, b00 = b0, c00 = c0;
        resampleStream.point = linePoint;
      }

      function ringEnd() {
        resampleLineTo(x0, y0, lambda0, a0, b0, c0, x00, y00, lambda00, a00, b00, c00, maxDepth, stream);
        resampleStream.lineEnd = lineEnd;
        lineEnd();
      }

      return resampleStream;
    };
  }

  var transformRadians = transformer({
    point: function(x, y) {
      this.stream.point(x * radians, y * radians);
    }
  });

  function transformRotate(rotate) {
    return transformer({
      point: function(x, y) {
        var r = rotate(x, y);
        return this.stream.point(r[0], r[1]);
      }
    });
  }

  function scaleTranslate(k, dx, dy) {
    function transform$$1(x, y) {
      return [dx + k * x, dy - k * y];
    }
    transform$$1.invert = function(x, y) {
      return [(x - dx) / k, (dy - y) / k];
    };
    return transform$$1;
  }

  function scaleTranslateRotate(k, dx, dy, alpha) {
    var cosAlpha = cos(alpha),
        sinAlpha = sin(alpha),
        a = cosAlpha * k,
        b = sinAlpha * k,
        ai = cosAlpha / k,
        bi = sinAlpha / k,
        ci = (sinAlpha * dy - cosAlpha * dx) / k,
        fi = (sinAlpha * dx + cosAlpha * dy) / k;
    function transform$$1(x, y) {
      return [a * x - b * y + dx, dy - b * x - a * y];
    }
    transform$$1.invert = function(x, y) {
      return [ai * x - bi * y + ci, fi - bi * x - ai * y];
    };
    return transform$$1;
  }

  function projection(project) {
    return projectionMutator(function() { return project; })();
  }

  function projectionMutator(projectAt) {
    var project,
        k = 150, // scale
        x = 480, y = 250, // translate
        lambda = 0, phi = 0, // center
        deltaLambda = 0, deltaPhi = 0, deltaGamma = 0, rotate, // pre-rotate
        alpha = 0, // post-rotate
        theta = null, preclip = clipAntimeridian, // pre-clip angle
        x0 = null, y0, x1, y1, postclip = identity$1, // post-clip extent
        delta2 = 0.5, // precision
        projectResample,
        projectTransform,
        projectRotateTransform,
        cache,
        cacheStream;

    function projection(point) {
      return projectRotateTransform(point[0] * radians, point[1] * radians);
    }

    function invert(point) {
      point = projectRotateTransform.invert(point[0], point[1]);
      return point && [point[0] * degrees, point[1] * degrees];
    }

    projection.stream = function(stream) {
      return cache && cacheStream === stream ? cache : cache = transformRadians(transformRotate(rotate)(preclip(projectResample(postclip(cacheStream = stream)))));
    };

    projection.preclip = function(_) {
      return arguments.length ? (preclip = _, theta = undefined, reset()) : preclip;
    };

    projection.postclip = function(_) {
      return arguments.length ? (postclip = _, x0 = y0 = x1 = y1 = null, reset()) : postclip;
    };

    projection.clipAngle = function(_) {
      return arguments.length ? (preclip = +_ ? clipCircle(theta = _ * radians) : (theta = null, clipAntimeridian), reset()) : theta * degrees;
    };

    projection.clipExtent = function(_) {
      return arguments.length ? (postclip = _ == null ? (x0 = y0 = x1 = y1 = null, identity$1) : clipRectangle(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
    };

    projection.scale = function(_) {
      return arguments.length ? (k = +_, recenter()) : k;
    };

    projection.translate = function(_) {
      return arguments.length ? (x = +_[0], y = +_[1], recenter()) : [x, y];
    };

    projection.center = function(_) {
      return arguments.length ? (lambda = _[0] % 360 * radians, phi = _[1] % 360 * radians, recenter()) : [lambda * degrees, phi * degrees];
    };

    projection.rotate = function(_) {
      return arguments.length ? (deltaLambda = _[0] % 360 * radians, deltaPhi = _[1] % 360 * radians, deltaGamma = _.length > 2 ? _[2] % 360 * radians : 0, recenter()) : [deltaLambda * degrees, deltaPhi * degrees, deltaGamma * degrees];
    };

    projection.angle = function(_) {
      return arguments.length ? (alpha = _ % 360 * radians, recenter()) : alpha * degrees;
    };

    projection.precision = function(_) {
      return arguments.length ? (projectResample = resample(projectTransform, delta2 = _ * _), reset()) : sqrt(delta2);
    };

    projection.fitExtent = function(extent, object) {
      return fitExtent(projection, extent, object);
    };

    projection.fitSize = function(size, object) {
      return fitSize(projection, size, object);
    };

    projection.fitWidth = function(width, object) {
      return fitWidth(projection, width, object);
    };

    projection.fitHeight = function(height, object) {
      return fitHeight(projection, height, object);
    };

    function recenter() {
      var center = scaleTranslateRotate(k, 0, 0, alpha).apply(null, project(lambda, phi)),
          transform$$1 = (alpha ? scaleTranslateRotate : scaleTranslate)(k, x - center[0], y - center[1], alpha);
      rotate = rotateRadians(deltaLambda, deltaPhi, deltaGamma);
      projectTransform = compose(project, transform$$1);
      projectRotateTransform = compose(rotate, projectTransform);
      projectResample = resample(projectTransform, delta2);
      return reset();
    }

    function reset() {
      cache = cacheStream = null;
      return projection;
    }

    return function() {
      project = projectAt.apply(this, arguments);
      projection.invert = project.invert && invert;
      return recenter();
    };
  }

  function azimuthalRaw(scale) {
    return function(x, y) {
      var cx = cos(x),
          cy = cos(y),
          k = scale(cx * cy);
      return [
        k * cy * sin(x),
        k * sin(y)
      ];
    }
  }

  function azimuthalInvert(angle) {
    return function(x, y) {
      var z = sqrt(x * x + y * y),
          c = angle(z),
          sc = sin(c),
          cc = cos(c);
      return [
        atan2(x * sc, z * cc),
        asin(z && y * sc / z)
      ];
    }
  }

  var azimuthalEqualAreaRaw = azimuthalRaw(function(cxcy) {
    return sqrt(2 / (1 + cxcy));
  });

  azimuthalEqualAreaRaw.invert = azimuthalInvert(function(z) {
    return 2 * asin(z / 2);
  });

  var azimuthalEquidistantRaw = azimuthalRaw(function(c) {
    return (c = acos(c)) && c / sin(c);
  });

  azimuthalEquidistantRaw.invert = azimuthalInvert(function(z) {
    return z;
  });

  function equirectangularRaw(lambda, phi) {
    return [lambda, phi];
  }

  equirectangularRaw.invert = equirectangularRaw;

  function geoEquirectangular() {
    return projection(equirectangularRaw)
        .scale(152.63);
  }

  function gnomonicRaw(x, y) {
    var cy = cos(y), k = cos(x) * cy;
    return [cy * sin(x) / k, sin(y) / k];
  }

  gnomonicRaw.invert = azimuthalInvert(atan);

  function gnomonic() {
    return projection(gnomonicRaw)
        .scale(144.049)
        .clipAngle(60);
  }

  function naturalEarth1Raw(lambda, phi) {
    var phi2 = phi * phi, phi4 = phi2 * phi2;
    return [
      lambda * (0.8707 - 0.131979 * phi2 + phi4 * (-0.013791 + phi4 * (0.003971 * phi2 - 0.001529 * phi4))),
      phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4)))
    ];
  }

  naturalEarth1Raw.invert = function(x, y) {
    var phi = y, i = 25, delta;
    do {
      var phi2 = phi * phi, phi4 = phi2 * phi2;
      phi -= delta = (phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4))) - y) /
          (1.007226 + phi2 * (0.015085 * 3 + phi4 * (-0.044475 * 7 + 0.028874 * 9 * phi2 - 0.005916 * 11 * phi4)));
    } while (abs(delta) > epsilon && --i > 0);
    return [
      x / (0.8707 + (phi2 = phi * phi) * (-0.131979 + phi2 * (-0.013791 + phi2 * phi2 * phi2 * (0.003971 - 0.001529 * phi2)))),
      phi
    ];
  };

  function naturalEarth1() {
    return projection(naturalEarth1Raw)
        .scale(175.295);
  }

  function orthographicRaw(x, y) {
    return [cos(y) * sin(x), sin(y)];
  }

  orthographicRaw.invert = azimuthalInvert(asin);

  function geoOrthographic() {
    return projection(orthographicRaw)
        .scale(249.5)
        .clipAngle(90 + epsilon);
  }

  var abs$1 = Math.abs;
  var atan$1 = Math.atan;
  var atan2$1 = Math.atan2;
  var cos$1 = Math.cos;
  var exp$1 = Math.exp;
  var floor$1 = Math.floor;
  var log$1 = Math.log;
  var max$1 = Math.max;
  var min$1 = Math.min;
  var pow$1 = Math.pow;
  var round = Math.round;
  var sign$1 = Math.sign || function(x) { return x > 0 ? 1 : x < 0 ? -1 : 0; };
  var sin$1 = Math.sin;
  var tan$1 = Math.tan;

  var epsilon$1 = 1e-6;
  var epsilon2$1 = 1e-12;
  var pi$1 = Math.PI;
  var halfPi$1 = pi$1 / 2;
  var quarterPi$1 = pi$1 / 4;
  var sqrt1_2 = Math.SQRT1_2;
  var sqrt2 = sqrt$1(2);
  var sqrtPi = sqrt$1(pi$1);
  var tau$1 = pi$1 * 2;
  var degrees$1 = 180 / pi$1;
  var radians$1 = pi$1 / 180;

  function sinci(x) {
    return x ? x / Math.sin(x) : 1;
  }

  function asin$1(x) {
    return x > 1 ? halfPi$1 : x < -1 ? -halfPi$1 : Math.asin(x);
  }

  function acos$1(x) {
    return x > 1 ? 0 : x < -1 ? pi$1 : Math.acos(x);
  }

  function sqrt$1(x) {
    return x > 0 ? Math.sqrt(x) : 0;
  }

  function tanh(x) {
    x = exp$1(2 * x);
    return (x - 1) / (x + 1);
  }

  function sinh(x) {
    return (exp$1(x) - exp$1(-x)) / 2;
  }

  function cosh(x) {
    return (exp$1(x) + exp$1(-x)) / 2;
  }

  function arsinh(x) {
    return log$1(x + sqrt$1(x * x + 1));
  }

  function arcosh(x) {
    return log$1(x + sqrt$1(x * x - 1));
  }

  function airyRaw(beta) {
    var tanBeta_2 = tan$1(beta / 2),
        b = 2 * log$1(cos$1(beta / 2)) / (tanBeta_2 * tanBeta_2);

    function forward(x, y) {
      var cosx = cos$1(x),
          cosy = cos$1(y),
          siny = sin$1(y),
          cosz = cosy * cosx,
          k = -((1 - cosz ? log$1((1 + cosz) / 2) / (1 - cosz) : -0.5) + b / (1 + cosz));
      return [k * cosy * sin$1(x), k * siny];
    }

    forward.invert = function(x, y) {
      var r = sqrt$1(x * x + y * y),
          z = -beta / 2,
          i = 50, delta;
      if (!r) return [0, 0];
      do {
        var z_2 = z / 2,
            cosz_2 = cos$1(z_2),
            sinz_2 = sin$1(z_2),
            tanz_2 = tan$1(z_2),
            lnsecz_2 = log$1(1 / cosz_2);
        z -= delta = (2 / tanz_2 * lnsecz_2 - b * tanz_2 - r) / (-lnsecz_2 / (sinz_2 * sinz_2) + 1 - b / (2 * cosz_2 * cosz_2));
      } while (abs$1(delta) > epsilon$1 && --i > 0);
      var sinz = sin$1(z);
      return [atan2$1(x * sinz, r * cos$1(z)), asin$1(y * sinz / r)];
    };

    return forward;
  }

  function airy() {
    var beta = halfPi$1,
        m = projectionMutator(airyRaw),
        p = m(beta);

    p.radius = function(_) {
      return arguments.length ? m(beta = _ * radians$1) : beta * degrees$1;
    };

    return p
        .scale(179.976)
        .clipAngle(147);
  }

  function aitoffRaw(x, y) {
    var cosy = cos$1(y), sincia = sinci(acos$1(cosy * cos$1(x /= 2)));
    return [2 * cosy * sin$1(x) * sincia, sin$1(y) * sincia];
  }

  // Abort if [x, y] is not within an ellipse centered at [0, 0] with
  // semi-major axis pi and semi-minor axis pi/2.
  aitoffRaw.invert = function(x, y) {
    if (x * x + 4 * y * y > pi$1 * pi$1 + epsilon$1) return;
    var x1 = x, y1 = y, i = 25;
    do {
      var sinx = sin$1(x1),
          sinx_2 = sin$1(x1 / 2),
          cosx_2 = cos$1(x1 / 2),
          siny = sin$1(y1),
          cosy = cos$1(y1),
          sin_2y = sin$1(2 * y1),
          sin2y = siny * siny,
          cos2y = cosy * cosy,
          sin2x_2 = sinx_2 * sinx_2,
          c = 1 - cos2y * cosx_2 * cosx_2,
          e = c ? acos$1(cosy * cosx_2) * sqrt$1(f = 1 / c) : f = 0,
          f,
          fx = 2 * e * cosy * sinx_2 - x,
          fy = e * siny - y,
          dxdx = f * (cos2y * sin2x_2 + e * cosy * cosx_2 * sin2y),
          dxdy = f * (0.5 * sinx * sin_2y - e * 2 * siny * sinx_2),
          dydx = f * 0.25 * (sin_2y * sinx_2 - e * siny * cos2y * sinx),
          dydy = f * (sin2y * cosx_2 + e * sin2x_2 * cosy),
          z = dxdy * dydx - dydy * dxdx;
      if (!z) break;
      var dx = (fy * dxdy - fx * dydy) / z,
          dy = (fx * dydx - fy * dxdx) / z;
      x1 -= dx, y1 -= dy;
    } while ((abs$1(dx) > epsilon$1 || abs$1(dy) > epsilon$1) && --i > 0);
    return [x1, y1];
  };

  function aitoff() {
    return projection(aitoffRaw)
        .scale(152.63);
  }

  function armadilloRaw(phi0) {
    var sinPhi0 = sin$1(phi0),
        cosPhi0 = cos$1(phi0),
        sPhi0 = phi0 >= 0 ? 1 : -1,
        tanPhi0 = tan$1(sPhi0 * phi0),
        k = (1 + sinPhi0 - cosPhi0) / 2;

    function forward(lambda, phi) {
      var cosPhi = cos$1(phi),
          cosLambda = cos$1(lambda /= 2);
      return [
        (1 + cosPhi) * sin$1(lambda),
        (sPhi0 * phi > -atan2$1(cosLambda, tanPhi0) - 1e-3 ? 0 : -sPhi0 * 10) + k + sin$1(phi) * cosPhi0 - (1 + cosPhi) * sinPhi0 * cosLambda // TODO D3 core should allow null or [NaN, NaN] to be returned.
      ];
    }

    forward.invert = function(x, y) {
      var lambda = 0,
          phi = 0,
          i = 50;
      do {
        var cosLambda = cos$1(lambda),
            sinLambda = sin$1(lambda),
            cosPhi = cos$1(phi),
            sinPhi = sin$1(phi),
            A = 1 + cosPhi,
            fx = A * sinLambda - x,
            fy = k + sinPhi * cosPhi0 - A * sinPhi0 * cosLambda - y,
            dxdLambda = A * cosLambda / 2,
            dxdPhi = -sinLambda * sinPhi,
            dydLambda = sinPhi0 * A * sinLambda / 2,
            dydPhi = cosPhi0 * cosPhi + sinPhi0 * cosLambda * sinPhi,
            denominator = dxdPhi * dydLambda - dydPhi * dxdLambda,
            dLambda = (fy * dxdPhi - fx * dydPhi) / denominator / 2,
            dPhi = (fx * dydLambda - fy * dxdLambda) / denominator;
        lambda -= dLambda, phi -= dPhi;
      } while ((abs$1(dLambda) > epsilon$1 || abs$1(dPhi) > epsilon$1) && --i > 0);
      return sPhi0 * phi > -atan2$1(cos$1(lambda), tanPhi0) - 1e-3 ? [lambda * 2, phi] : null;
    };

    return forward;
  }

  function armadillo() {
    var phi0 = 20 * radians$1,
        sPhi0 = phi0 >= 0 ? 1 : -1,
        tanPhi0 = tan$1(sPhi0 * phi0),
        m = projectionMutator(armadilloRaw),
        p = m(phi0),
        stream_ = p.stream;

    p.parallel = function(_) {
      if (!arguments.length) return phi0 * degrees$1;
      tanPhi0 = tan$1((sPhi0 = (phi0 = _ * radians$1) >= 0 ? 1 : -1) * phi0);
      return m(phi0);
    };

    p.stream = function(stream) {
      var rotate = p.rotate(),
          rotateStream = stream_(stream),
          sphereStream = (p.rotate([0, 0]), stream_(stream));
      p.rotate(rotate);
      rotateStream.sphere = function() {
        sphereStream.polygonStart(), sphereStream.lineStart();
        for (var lambda = sPhi0 * -180; sPhi0 * lambda < 180; lambda += sPhi0 * 90) sphereStream.point(lambda, sPhi0 * 90);
        while (sPhi0 * (lambda -= phi0) >= -180) { // TODO precision?
          sphereStream.point(lambda, sPhi0 * -atan2$1(cos$1(lambda * radians$1 / 2), tanPhi0) * degrees$1);
        }
        sphereStream.lineEnd(), sphereStream.polygonEnd();
      };
      return rotateStream;
    };

    return p
        .scale(218.695)
        .center([0, 28.0974]);
  }

  function augustRaw(lambda, phi) {
    var tanPhi = tan$1(phi / 2),
        k = sqrt$1(1 - tanPhi * tanPhi),
        c = 1 + k * cos$1(lambda /= 2),
        x = sin$1(lambda) * k / c,
        y = tanPhi / c,
        x2 = x * x,
        y2 = y * y;
    return [
      4 / 3 * x * (3 + x2 - 3 * y2),
      4 / 3 * y * (3 + 3 * x2 - y2)
    ];
  }

  augustRaw.invert = function(x, y) {
    x *= 3 / 8, y *= 3 / 8;
    if (!x && abs$1(y) > 1) return null;
    var x2 = x * x,
        y2 = y * y,
        s = 1 + x2 + y2,
        sin3Eta = sqrt$1((s - sqrt$1(s * s - 4 * y * y)) / 2),
        eta = asin$1(sin3Eta) / 3,
        xi = sin3Eta ? arcosh(abs$1(y / sin3Eta)) / 3 : arsinh(abs$1(x)) / 3,
        cosEta = cos$1(eta),
        coshXi = cosh(xi),
        d = coshXi * coshXi - cosEta * cosEta;
    return [
      sign$1(x) * 2 * atan2$1(sinh(xi) * cosEta, 0.25 - d),
      sign$1(y) * 2 * atan2$1(coshXi * sin$1(eta), 0.25 + d)
    ];
  };

  function august() {
    return projection(augustRaw)
        .scale(66.1603);
  }

  var sqrt8 = sqrt$1(8),
      phi0$1 = log$1(1 + sqrt2);

  function bakerRaw(lambda, phi) {
    var phi0 = abs$1(phi);
    return phi0 < quarterPi$1
        ? [lambda, log$1(tan$1(quarterPi$1 + phi / 2))]
        : [lambda * cos$1(phi0) * (2 * sqrt2 - 1 / sin$1(phi0)), sign$1(phi) * (2 * sqrt2 * (phi0 - quarterPi$1) - log$1(tan$1(phi0 / 2)))];
  }

  bakerRaw.invert = function(x, y) {
    if ((y0 = abs$1(y)) < phi0$1) return [x, 2 * atan$1(exp$1(y)) - halfPi$1];
    var phi = quarterPi$1, i = 25, delta, y0;
    do {
      var cosPhi_2 = cos$1(phi / 2), tanPhi_2 = tan$1(phi / 2);
      phi -= delta = (sqrt8 * (phi - quarterPi$1) - log$1(tanPhi_2) - y0) / (sqrt8 - cosPhi_2 * cosPhi_2 / (2 * tanPhi_2));
    } while (abs$1(delta) > epsilon2$1 && --i > 0);
    return [x / (cos$1(phi) * (sqrt8 - 1 / sin$1(phi))), sign$1(y) * phi];
  };

  function baker() {
    return projection(bakerRaw)
        .scale(112.314);
  }

  function berghausRaw(lobes) {
    var k = 2 * pi$1 / lobes;

    function forward(lambda, phi) {
      var p = azimuthalEquidistantRaw(lambda, phi);
      if (abs$1(lambda) > halfPi$1) { // back hemisphere
        var theta = atan2$1(p[1], p[0]),
            r = sqrt$1(p[0] * p[0] + p[1] * p[1]),
            theta0 = k * round((theta - halfPi$1) / k) + halfPi$1,
            alpha = atan2$1(sin$1(theta -= theta0), 2 - cos$1(theta)); // angle relative to lobe end
        theta = theta0 + asin$1(pi$1 / r * sin$1(alpha)) - alpha;
        p[0] = r * cos$1(theta);
        p[1] = r * sin$1(theta);
      }
      return p;
    }

    forward.invert = function(x, y) {
      var r = sqrt$1(x * x + y * y);
      if (r > halfPi$1) {
        var theta = atan2$1(y, x),
            theta0 = k * round((theta - halfPi$1) / k) + halfPi$1,
            s = theta > theta0 ? -1 : 1,
            A = r * cos$1(theta0 - theta),
            cotAlpha = 1 / tan$1(s * acos$1((A - pi$1) / sqrt$1(pi$1 * (pi$1 - 2 * A) + r * r)));
        theta = theta0 + 2 * atan$1((cotAlpha + s * sqrt$1(cotAlpha * cotAlpha - 3)) / 3);
        x = r * cos$1(theta), y = r * sin$1(theta);
      }
      return azimuthalEquidistantRaw.invert(x, y);
    };

    return forward;
  }

  function berghaus() {
    var lobes = 5,
        m = projectionMutator(berghausRaw),
        p = m(lobes),
        projectionStream = p.stream,
        epsilon = 1e-2,
        cr = -cos$1(epsilon * radians$1),
        sr = sin$1(epsilon * radians$1);

    p.lobes = function(_) {
      return arguments.length ? m(lobes = +_) : lobes;
    };

    p.stream = function(stream) {
      var rotate = p.rotate(),
          rotateStream = projectionStream(stream),
          sphereStream = (p.rotate([0, 0]), projectionStream(stream));
      p.rotate(rotate);
      rotateStream.sphere = function() {
        sphereStream.polygonStart(), sphereStream.lineStart();
        for (var i = 0, delta = 360 / lobes, delta0 = 2 * pi$1 / lobes, phi = 90 - 180 / lobes, phi0 = halfPi$1; i < lobes; ++i, phi -= delta, phi0 -= delta0) {
          sphereStream.point(atan2$1(sr * cos$1(phi0), cr) * degrees$1, asin$1(sr * sin$1(phi0)) * degrees$1);
          if (phi < -90) {
            sphereStream.point(-90, -180 - phi - epsilon);
            sphereStream.point(-90, -180 - phi + epsilon);
          } else {
            sphereStream.point(90, phi + epsilon);
            sphereStream.point(90, phi - epsilon);
          }
        }
        sphereStream.lineEnd(), sphereStream.polygonEnd();
      };
      return rotateStream;
    };

    return p
        .scale(87.8076)
        .center([0, 17.1875])
        .clipAngle(180 - 1e-3);
  }

  function hammerRaw(A, B) {
    if (arguments.length < 2) B = A;
    if (B === 1) return azimuthalEqualAreaRaw;
    if (B === Infinity) return hammerQuarticAuthalicRaw;

    function forward(lambda, phi) {
      var coordinates = azimuthalEqualAreaRaw(lambda / B, phi);
      coordinates[0] *= A;
      return coordinates;
    }

    forward.invert = function(x, y) {
      var coordinates = azimuthalEqualAreaRaw.invert(x / A, y);
      coordinates[0] *= B;
      return coordinates;
    };

    return forward;
  }

  function hammerQuarticAuthalicRaw(lambda, phi) {
    return [
      lambda * cos$1(phi) / cos$1(phi /= 2),
      2 * sin$1(phi)
    ];
  }

  hammerQuarticAuthalicRaw.invert = function(x, y) {
    var phi = 2 * asin$1(y / 2);
    return [
      x * cos$1(phi / 2) / cos$1(phi),
      phi
    ];
  };

  function hammer() {
    var B = 2,
        m = projectionMutator(hammerRaw),
        p = m(B);

    p.coefficient = function(_) {
      if (!arguments.length) return B;
      return m(B = +_);
    };

    return p
      .scale(169.529);
  }

  // Bertin 1953 as a modified Briesemeister
  // https://bl.ocks.org/Fil/5b9ee9636dfb6ffa53443c9006beb642
  function bertin1953Raw() {
    var hammer$$1 = hammerRaw(1.68, 2),
        fu = 1.4, k = 12;

    return function(lambda, phi) {

      if (lambda + phi < -fu) {
        var u = (lambda - phi + 1.6) * (lambda + phi + fu) / 8;
        lambda += u;
        phi -= 0.8 * u * sin$1(phi + pi$1 / 2);
      }

      var r = hammer$$1(lambda, phi);

      var d = (1 - cos$1(lambda * phi)) / k;

      if (r[1] < 0) {
        r[0] *= 1 + d;
      }
      if (r[1] > 0) {
        r[1] *= 1 + d / 1.5 * r[0] * r[0];
      }

      return r;
    };
  }

  function bertin() {
    var p = projection(bertin1953Raw());

    p.rotate([-16.5, -42]);
    delete p.rotate;

    return p
      .scale(176.57)
      .center([7.93, 0.09]);
  }

  function mollweideBromleyTheta(cp, phi) {
    var cpsinPhi = cp * sin$1(phi), i = 30, delta;
    do phi -= delta = (phi + sin$1(phi) - cpsinPhi) / (1 + cos$1(phi));
    while (abs$1(delta) > epsilon$1 && --i > 0);
    return phi / 2;
  }

  function mollweideBromleyRaw(cx, cy, cp) {

    function forward(lambda, phi) {
      return [cx * lambda * cos$1(phi = mollweideBromleyTheta(cp, phi)), cy * sin$1(phi)];
    }

    forward.invert = function(x, y) {
      return y = asin$1(y / cy), [x / (cx * cos$1(y)), asin$1((2 * y + sin$1(2 * y)) / cp)];
    };

    return forward;
  }

  var mollweideRaw = mollweideBromleyRaw(sqrt2 / halfPi$1, sqrt2, pi$1);

  function mollweide() {
    return projection(mollweideRaw)
        .scale(169.529);
  }

  var k = 2.00276,
      w = 1.11072;

  function boggsRaw(lambda, phi) {
    var theta = mollweideBromleyTheta(pi$1, phi);
    return [k * lambda / (1 / cos$1(phi) + w / cos$1(theta)), (phi + sqrt2 * sin$1(theta)) / k];
  }

  boggsRaw.invert = function(x, y) {
    var ky = k * y, theta = y < 0 ? -quarterPi$1 : quarterPi$1, i = 25, delta, phi;
    do {
      phi = ky - sqrt2 * sin$1(theta);
      theta -= delta = (sin$1(2 * theta) + 2 * theta - pi$1 * sin$1(phi)) / (2 * cos$1(2 * theta) + 2 + pi$1 * cos$1(phi) * sqrt2 * cos$1(theta));
    } while (abs$1(delta) > epsilon$1 && --i > 0);
    phi = ky - sqrt2 * sin$1(theta);
    return [x * (1 / cos$1(phi) + w / cos$1(theta)) / k, phi];
  };

  function boggs() {
    return projection(boggsRaw)
        .scale(160.857);
  }

  function parallel1(projectAt) {
    var phi0 = 0,
        m = projectionMutator(projectAt),
        p = m(phi0);

    p.parallel = function(_) {
      return arguments.length ? m(phi0 = _ * radians$1) : phi0 * degrees$1;
    };

    return p;
  }

  function sinusoidalRaw(lambda, phi) {
    return [lambda * cos$1(phi), phi];
  }

  sinusoidalRaw.invert = function(x, y) {
    return [x / cos$1(y), y];
  };

  function sinusoidal() {
    return projection(sinusoidalRaw)
        .scale(152.63);
  }

  function bonneRaw(phi0) {
    if (!phi0) return sinusoidalRaw;
    var cotPhi0 = 1 / tan$1(phi0);

    function forward(lambda, phi) {
      var rho = cotPhi0 + phi0 - phi,
          e = rho ? lambda * cos$1(phi) / rho : rho;
      return [rho * sin$1(e), cotPhi0 - rho * cos$1(e)];
    }

    forward.invert = function(x, y) {
      var rho = sqrt$1(x * x + (y = cotPhi0 - y) * y),
          phi = cotPhi0 + phi0 - rho;
      return [rho / cos$1(phi) * atan2$1(x, y), phi];
    };

    return forward;
  }

  function bonne() {
    return parallel1(bonneRaw)
        .scale(123.082)
        .center([0, 26.1441])
        .parallel(45);
  }

  function bottomleyRaw(sinPsi) {

    function forward(lambda, phi) {
      var rho = halfPi$1 - phi,
          eta = rho ? lambda * sinPsi * sin$1(rho) / rho : rho;
      return [rho * sin$1(eta) / sinPsi, halfPi$1 - rho * cos$1(eta)];
    }

    forward.invert = function(x, y) {
      var x1 = x * sinPsi,
          y1 = halfPi$1 - y,
          rho = sqrt$1(x1 * x1 + y1 * y1),
          eta = atan2$1(x1, y1);
      return [(rho ? rho / sin$1(rho) : 1) * eta / sinPsi, halfPi$1 - rho];
    };

    return forward;
  }

  function bottomley() {
    var sinPsi = 0.5,
        m = projectionMutator(bottomleyRaw),
        p = m(sinPsi);

    p.fraction = function(_) {
      return arguments.length ? m(sinPsi = +_) : sinPsi;
    };

    return p
        .scale(158.837);
  }

  var bromleyRaw = mollweideBromleyRaw(1, 4 / pi$1, pi$1);

  function bromley() {
    return projection(bromleyRaw)
        .scale(152.63);
  }

  // Azimuthal distance.
  function distance$1(dPhi, c1, s1, c2, s2, dLambda) {
    var cosdLambda = cos$1(dLambda), r;
    if (abs$1(dPhi) > 1 || abs$1(dLambda) > 1) {
      r = acos$1(s1 * s2 + c1 * c2 * cosdLambda);
    } else {
      var sindPhi = sin$1(dPhi / 2), sindLambda = sin$1(dLambda / 2);
      r = 2 * asin$1(sqrt$1(sindPhi * sindPhi + c1 * c2 * sindLambda * sindLambda));
    }
    return abs$1(r) > epsilon$1 ? [r, atan2$1(c2 * sin$1(dLambda), c1 * s2 - s1 * c2 * cosdLambda)] : [0, 0];
  }

  // Angle opposite a, and contained between sides of lengths b and c.
  function angle$1(b, c, a) {
    return acos$1((b * b + c * c - a * a) / (2 * b * c));
  }

  // Normalize longitude.
  function longitude(lambda) {
    return lambda - 2 * pi$1 * floor$1((lambda + pi$1) / (2 * pi$1));
  }

  function chamberlinRaw(p0, p1, p2) {
    var points = [
      [p0[0], p0[1], sin$1(p0[1]), cos$1(p0[1])],
      [p1[0], p1[1], sin$1(p1[1]), cos$1(p1[1])],
      [p2[0], p2[1], sin$1(p2[1]), cos$1(p2[1])]
    ];

    for (var a = points[2], b, i = 0; i < 3; ++i, a = b) {
      b = points[i];
      a.v = distance$1(b[1] - a[1], a[3], a[2], b[3], b[2], b[0] - a[0]);
      a.point = [0, 0];
    }

    var beta0 = angle$1(points[0].v[0], points[2].v[0], points[1].v[0]),
        beta1 = angle$1(points[0].v[0], points[1].v[0], points[2].v[0]),
        beta2 = pi$1 - beta0;

    points[2].point[1] = 0;
    points[0].point[0] = -(points[1].point[0] = points[0].v[0] / 2);

    var mean = [
      points[2].point[0] = points[0].point[0] + points[2].v[0] * cos$1(beta0),
      2 * (points[0].point[1] = points[1].point[1] = points[2].v[0] * sin$1(beta0))
    ];

    function forward(lambda, phi) {
      var sinPhi = sin$1(phi),
          cosPhi = cos$1(phi),
          v = new Array(3), i;

      // Compute distance and azimuth from control points.
      for (i = 0; i < 3; ++i) {
        var p = points[i];
        v[i] = distance$1(phi - p[1], p[3], p[2], cosPhi, sinPhi, lambda - p[0]);
        if (!v[i][0]) return p.point;
        v[i][1] = longitude(v[i][1] - p.v[1]);
      }

      // Arithmetic mean of interception points.
      var point = mean.slice();
      for (i = 0; i < 3; ++i) {
        var j = i == 2 ? 0 : i + 1;
        var a = angle$1(points[i].v[0], v[i][0], v[j][0]);
        if (v[i][1] < 0) a = -a;

        if (!i) {
          point[0] += v[i][0] * cos$1(a);
          point[1] -= v[i][0] * sin$1(a);
        } else if (i == 1) {
          a = beta1 - a;
          point[0] -= v[i][0] * cos$1(a);
          point[1] -= v[i][0] * sin$1(a);
        } else {
          a = beta2 - a;
          point[0] += v[i][0] * cos$1(a);
          point[1] += v[i][0] * sin$1(a);
        }
      }

      point[0] /= 3, point[1] /= 3;
      return point;
    }

    return forward;
  }

  function pointRadians$1(p) {
    return p[0] *= radians$1, p[1] *= radians$1, p;
  }

  function chamberlinAfrica() {
    return chamberlin([0, 22], [45, 22], [22.5, -22])
        .scale(380)
        .center([22.5, 2]);
  }

  function chamberlin(p0, p1, p2) { // TODO order matters!
    var c = centroid({type: "MultiPoint", coordinates: [p0, p1, p2]}),
        R = [-c[0], -c[1]],
        r = rotation(R),
        p = projection(chamberlinRaw(pointRadians$1(r(p0)), pointRadians$1(r(p1)), pointRadians$1(r(p2)))).rotate(R),
        center = p.center;

    delete p.rotate;

    p.center = function(_) {
      return arguments.length ? center(r(_)) : r.invert(center());
    };

    return p
        .clipAngle(90);
  }

  function collignonRaw(lambda, phi) {
    var alpha = sqrt$1(1 - sin$1(phi));
    return [(2 / sqrtPi) * lambda * alpha, sqrtPi * (1 - alpha)];
  }

  collignonRaw.invert = function(x, y) {
    var lambda = (lambda = y / sqrtPi - 1) * lambda;
    return [lambda > 0 ? x * sqrt$1(pi$1 / lambda) / 2 : 0, asin$1(1 - lambda)];
  };

  function collignon() {
    return projection(collignonRaw)
        .scale(95.6464)
        .center([0, 30]);
  }

  function craigRaw(phi0) {
    var tanPhi0 = tan$1(phi0);

    function forward(lambda, phi) {
      return [lambda, (lambda ? lambda / sin$1(lambda) : 1) * (sin$1(phi) * cos$1(lambda) - tanPhi0 * cos$1(phi))];
    }

    forward.invert = tanPhi0 ? function(x, y) {
      if (x) y *= sin$1(x) / x;
      var cosLambda = cos$1(x);
      return [x, 2 * atan2$1(sqrt$1(cosLambda * cosLambda + tanPhi0 * tanPhi0 - y * y) - cosLambda, tanPhi0 - y)];
    } : function(x, y) {
      return [x, asin$1(x ? y * tan$1(x) / x : y)];
    };

    return forward;
  }

  function craig() {
    return parallel1(craigRaw)
        .scale(249.828)
        .clipAngle(90);
  }

  var sqrt3 = sqrt$1(3);

  function crasterRaw(lambda, phi) {
    return [sqrt3 * lambda * (2 * cos$1(2 * phi / 3) - 1) / sqrtPi, sqrt3 * sqrtPi * sin$1(phi / 3)];
  }

  crasterRaw.invert = function(x, y) {
    var phi = 3 * asin$1(y / (sqrt3 * sqrtPi));
    return [sqrtPi * x / (sqrt3 * (2 * cos$1(2 * phi / 3) - 1)), phi];
  };

  function craster() {
    return projection(crasterRaw)
        .scale(156.19);
  }

  function cylindricalEqualAreaRaw$1(phi0) {
    var cosPhi0 = cos$1(phi0);

    function forward(lambda, phi) {
      return [lambda * cosPhi0, sin$1(phi) / cosPhi0];
    }

    forward.invert = function(x, y) {
      return [x / cosPhi0, asin$1(y * cosPhi0)];
    };

    return forward;
  }

  function cylindricalEqualArea() {
    return parallel1(cylindricalEqualAreaRaw$1)
        .parallel(38.58) // acos(sqrt(width / height / pi)) * radians
        .scale(195.044); // width / (sqrt(width / height / pi) * 2 * pi)
  }

  function cylindricalStereographicRaw(phi0) {
    var cosPhi0 = cos$1(phi0);

    function forward(lambda, phi) {
      return [lambda * cosPhi0, (1 + cosPhi0) * tan$1(phi / 2)];
    }

    forward.invert = function(x, y) {
      return [x / cosPhi0, atan$1(y / (1 + cosPhi0)) * 2];
    };

    return forward;
  }

  function cylindricalStereographic() {
    return parallel1(cylindricalStereographicRaw)
        .scale(124.75);
  }

  function eckert1Raw(lambda, phi) {
    var alpha = sqrt$1(8 / (3 * pi$1));
    return [
      alpha * lambda * (1 - abs$1(phi) / pi$1),
      alpha * phi
    ];
  }

  eckert1Raw.invert = function(x, y) {
    var alpha = sqrt$1(8 / (3 * pi$1)),
        phi = y / alpha;
    return [
      x / (alpha * (1 - abs$1(phi) / pi$1)),
      phi
    ];
  };

  function eckert1() {
    return projection(eckert1Raw)
        .scale(165.664);
  }

  function eckert2Raw(lambda, phi) {
    var alpha = sqrt$1(4 - 3 * sin$1(abs$1(phi)));
    return [
      2 / sqrt$1(6 * pi$1) * lambda * alpha,
      sign$1(phi) * sqrt$1(2 * pi$1 / 3) * (2 - alpha)
    ];
  }

  eckert2Raw.invert = function(x, y) {
    var alpha = 2 - abs$1(y) / sqrt$1(2 * pi$1 / 3);
    return [
      x * sqrt$1(6 * pi$1) / (2 * alpha),
      sign$1(y) * asin$1((4 - alpha * alpha) / 3)
    ];
  };

  function eckert2() {
    return projection(eckert2Raw)
        .scale(165.664);
  }

  function eckert3Raw(lambda, phi) {
    var k = sqrt$1(pi$1 * (4 + pi$1));
    return [
      2 / k * lambda * (1 + sqrt$1(1 - 4 * phi * phi / (pi$1 * pi$1))),
      4 / k * phi
    ];
  }

  eckert3Raw.invert = function(x, y) {
    var k = sqrt$1(pi$1 * (4 + pi$1)) / 2;
    return [
      x * k / (1 + sqrt$1(1 - y * y * (4 + pi$1) / (4 * pi$1))),
      y * k / 2
    ];
  };

  function eckert3() {
    return projection(eckert3Raw)
        .scale(180.739);
  }

  function eckert4Raw(lambda, phi) {
    var k = (2 + halfPi$1) * sin$1(phi);
    phi /= 2;
    for (var i = 0, delta = Infinity; i < 10 && abs$1(delta) > epsilon$1; i++) {
      var cosPhi = cos$1(phi);
      phi -= delta = (phi + sin$1(phi) * (cosPhi + 2) - k) / (2 * cosPhi * (1 + cosPhi));
    }
    return [
      2 / sqrt$1(pi$1 * (4 + pi$1)) * lambda * (1 + cos$1(phi)),
      2 * sqrt$1(pi$1 / (4 + pi$1)) * sin$1(phi)
    ];
  }

  eckert4Raw.invert = function(x, y) {
    var A = y * sqrt$1((4 + pi$1) / pi$1) / 2,
        k = asin$1(A),
        c = cos$1(k);
    return [
      x / (2 / sqrt$1(pi$1 * (4 + pi$1)) * (1 + c)),
      asin$1((k + A * (c + 2)) / (2 + halfPi$1))
    ];
  };

  function eckert4() {
    return projection(eckert4Raw)
        .scale(180.739);
  }

  function eckert5Raw(lambda, phi) {
    return [
      lambda * (1 + cos$1(phi)) / sqrt$1(2 + pi$1),
      2 * phi / sqrt$1(2 + pi$1)
    ];
  }

  eckert5Raw.invert = function(x, y) {
    var k = sqrt$1(2 + pi$1),
        phi = y * k / 2;
    return [
      k * x / (1 + cos$1(phi)),
      phi
    ];
  };

  function eckert5() {
    return projection(eckert5Raw)
        .scale(173.044);
  }

  function eckert6Raw(lambda, phi) {
    var k = (1 + halfPi$1) * sin$1(phi);
    for (var i = 0, delta = Infinity; i < 10 && abs$1(delta) > epsilon$1; i++) {
      phi -= delta = (phi + sin$1(phi) - k) / (1 + cos$1(phi));
    }
    k = sqrt$1(2 + pi$1);
    return [
      lambda * (1 + cos$1(phi)) / k,
      2 * phi / k
    ];
  }

  eckert6Raw.invert = function(x, y) {
    var j = 1 + halfPi$1,
        k = sqrt$1(j / 2);
    return [
      x * 2 * k / (1 + cos$1(y *= k)),
      asin$1((y + sin$1(y)) / j)
    ];
  };

  function eckert6() {
    return projection(eckert6Raw)
        .scale(173.044);
  }

  var eisenlohrK = 3 + 2 * sqrt2;

  function eisenlohrRaw(lambda, phi) {
    var s0 = sin$1(lambda /= 2),
        c0 = cos$1(lambda),
        k = sqrt$1(cos$1(phi)),
        c1 = cos$1(phi /= 2),
        t = sin$1(phi) / (c1 + sqrt2 * c0 * k),
        c = sqrt$1(2 / (1 + t * t)),
        v = sqrt$1((sqrt2 * c1 + (c0 + s0) * k) / (sqrt2 * c1 + (c0 - s0) * k));
    return [
      eisenlohrK * (c * (v - 1 / v) - 2 * log$1(v)),
      eisenlohrK * (c * t * (v + 1 / v) - 2 * atan$1(t))
    ];
  }

  eisenlohrRaw.invert = function(x, y) {
    if (!(p = augustRaw.invert(x / 1.2, y * 1.065))) return null;
    var lambda = p[0], phi = p[1], i = 20, p;
    x /= eisenlohrK, y /= eisenlohrK;
    do {
      var _0 = lambda / 2,
          _1 = phi / 2,
          s0 = sin$1(_0),
          c0 = cos$1(_0),
          s1 = sin$1(_1),
          c1 = cos$1(_1),
          cos1 = cos$1(phi),
          k = sqrt$1(cos1),
          t = s1 / (c1 + sqrt2 * c0 * k),
          t2 = t * t,
          c = sqrt$1(2 / (1 + t2)),
          v0 = (sqrt2 * c1 + (c0 + s0) * k),
          v1 = (sqrt2 * c1 + (c0 - s0) * k),
          v2 = v0 / v1,
          v = sqrt$1(v2),
          vm1v = v - 1 / v,
          vp1v = v + 1 / v,
          fx = c * vm1v - 2 * log$1(v) - x,
          fy = c * t * vp1v - 2 * atan$1(t) - y,
          deltatDeltaLambda = s1 && sqrt1_2 * k * s0 * t2 / s1,
          deltatDeltaPhi = (sqrt2 * c0 * c1 + k) / (2 * (c1 + sqrt2 * c0 * k) * (c1 + sqrt2 * c0 * k) * k),
          deltacDeltat = -0.5 * t * c * c * c,
          deltacDeltaLambda = deltacDeltat * deltatDeltaLambda,
          deltacDeltaPhi = deltacDeltat * deltatDeltaPhi,
          A = (A = 2 * c1 + sqrt2 * k * (c0 - s0)) * A * v,
          deltavDeltaLambda = (sqrt2 * c0 * c1 * k + cos1) / A,
          deltavDeltaPhi = -(sqrt2 * s0 * s1) / (k * A),
          deltaxDeltaLambda = vm1v * deltacDeltaLambda - 2 * deltavDeltaLambda / v + c * (deltavDeltaLambda + deltavDeltaLambda / v2),
          deltaxDeltaPhi = vm1v * deltacDeltaPhi - 2 * deltavDeltaPhi / v + c * (deltavDeltaPhi + deltavDeltaPhi / v2),
          deltayDeltaLambda = t * vp1v * deltacDeltaLambda - 2 * deltatDeltaLambda / (1 + t2) + c * vp1v * deltatDeltaLambda + c * t * (deltavDeltaLambda - deltavDeltaLambda / v2),
          deltayDeltaPhi = t * vp1v * deltacDeltaPhi - 2 * deltatDeltaPhi / (1 + t2) + c * vp1v * deltatDeltaPhi + c * t * (deltavDeltaPhi - deltavDeltaPhi / v2),
          denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda;
      if (!denominator) break;
      var deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator,
          deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
      lambda -= deltaLambda;
      phi = max$1(-halfPi$1, min$1(halfPi$1, phi - deltaPhi));
    } while ((abs$1(deltaLambda) > epsilon$1 || abs$1(deltaPhi) > epsilon$1) && --i > 0);
    return abs$1(abs$1(phi) - halfPi$1) < epsilon$1 ? [0, phi] : i && [lambda, phi];
  };

  function eisenlohr() {
    return projection(eisenlohrRaw)
        .scale(62.5271);
  }

  var faheyK = cos$1(35 * radians$1);

  function faheyRaw(lambda, phi) {
    var t = tan$1(phi / 2);
    return [lambda * faheyK * sqrt$1(1 - t * t), (1 + faheyK) * t];
  }

  faheyRaw.invert = function(x, y) {
    var t = y / (1 + faheyK);
    return [x && x / (faheyK * sqrt$1(1 - t * t)), 2 * atan$1(t)];
  };

  function fahey() {
    return projection(faheyRaw)
        .scale(137.152);
  }

  function foucautRaw(lambda, phi) {
    var k = phi / 2, cosk = cos$1(k);
    return [ 2 * lambda / sqrtPi * cos$1(phi) * cosk * cosk, sqrtPi * tan$1(k)];
  }

  foucautRaw.invert = function(x, y) {
    var k = atan$1(y / sqrtPi), cosk = cos$1(k), phi = 2 * k;
    return [x * sqrtPi / 2 / (cos$1(phi) * cosk * cosk), phi];
  };

  function foucaut() {
    return projection(foucautRaw)
        .scale(135.264);
  }

  function gilbertForward(point) {
    return [point[0] / 2, asin$1(tan$1(point[1] / 2 * radians$1)) * degrees$1];
  }

  function gilbertInvert(point) {
    return [point[0] * 2, 2 * atan$1(sin$1(point[1] * radians$1)) * degrees$1];
  }

  function gilbert(projectionType) {
    if (projectionType == null) projectionType = geoOrthographic;
    var projection$$1 = projectionType(),
        equirectangular = geoEquirectangular().scale(degrees$1).precision(0).clipAngle(null).translate([0, 0]); // antimeridian cutting

    function gilbert(point) {
      return projection$$1(gilbertForward(point));
    }

    if (projection$$1.invert) gilbert.invert = function(point) {
      return gilbertInvert(projection$$1.invert(point));
    };

    gilbert.stream = function(stream) {
      var s1 = projection$$1.stream(stream), s0 = equirectangular.stream({
        point: function(lambda, phi) { s1.point(lambda / 2, asin$1(tan$1(-phi / 2 * radians$1)) * degrees$1); },
        lineStart: function() { s1.lineStart(); },
        lineEnd: function() { s1.lineEnd(); },
        polygonStart: function() { s1.polygonStart(); },
        polygonEnd: function() { s1.polygonEnd(); }
      });
      s0.sphere = s1.sphere;
      return s0;
    };

    function property(name) {
      gilbert[name] = function(_) {
        return arguments.length ? (projection$$1[name](_), gilbert) : projection$$1[name]();
      };
    }

    gilbert.rotate = function(_) {
      return arguments.length ? (equirectangular.rotate(_), gilbert) : equirectangular.rotate();
    };

    gilbert.center = function(_) {
      return arguments.length ? (projection$$1.center(gilbertForward(_)), gilbert) : gilbertInvert(projection$$1.center());
    };

    property("angle");
    property("clipAngle");
    property("clipExtent");
    property("scale");
    property("translate");
    property("precision");

    return gilbert
        .scale(249.5);
  }

  function gingeryRaw(rho, n) {
    var k = 2 * pi$1 / n,
        rho2 = rho * rho;

    function forward(lambda, phi) {
      var p = azimuthalEquidistantRaw(lambda, phi),
          x = p[0],
          y = p[1],
          r2 = x * x + y * y;

      if (r2 > rho2) {
        var r = sqrt$1(r2),
            theta = atan2$1(y, x),
            theta0 = k * round(theta / k),
            alpha = theta - theta0,
            rhoCosAlpha = rho * cos$1(alpha),
            k_ = (rho * sin$1(alpha) - alpha * sin$1(rhoCosAlpha)) / (halfPi$1 - rhoCosAlpha),
            s_ = gingeryLength(alpha, k_),
            e = (pi$1 - rho) / gingeryIntegrate(s_, rhoCosAlpha, pi$1);

        x = r;
        var i = 50, delta;
        do {
          x -= delta = (rho + gingeryIntegrate(s_, rhoCosAlpha, x) * e - r) / (s_(x) * e);
        } while (abs$1(delta) > epsilon$1 && --i > 0);

        y = alpha * sin$1(x);
        if (x < halfPi$1) y -= k_ * (x - halfPi$1);

        var s = sin$1(theta0),
            c = cos$1(theta0);
        p[0] = x * c - y * s;
        p[1] = x * s + y * c;
      }
      return p;
    }

    forward.invert = function(x, y) {
      var r2 = x * x + y * y;
      if (r2 > rho2) {
        var r = sqrt$1(r2),
            theta = atan2$1(y, x),
            theta0 = k * round(theta / k),
            dTheta = theta - theta0;

        x = r * cos$1(dTheta);
        y = r * sin$1(dTheta);

        var x_halfPi = x - halfPi$1,
            sinx = sin$1(x),
            alpha = y / sinx,
            delta = x < halfPi$1 ? Infinity : 0,
            i = 10;

        while (true) {
          var rhosinAlpha = rho * sin$1(alpha),
              rhoCosAlpha = rho * cos$1(alpha),
              sinRhoCosAlpha = sin$1(rhoCosAlpha),
              halfPi_RhoCosAlpha = halfPi$1 - rhoCosAlpha,
              k_ = (rhosinAlpha - alpha * sinRhoCosAlpha) / halfPi_RhoCosAlpha,
              s_ = gingeryLength(alpha, k_);

          if (abs$1(delta) < epsilon2$1 || !--i) break;

          alpha -= delta = (alpha * sinx - k_ * x_halfPi - y) / (
            sinx - x_halfPi * 2 * (
              halfPi_RhoCosAlpha * (rhoCosAlpha + alpha * rhosinAlpha * cos$1(rhoCosAlpha) - sinRhoCosAlpha) -
              rhosinAlpha * (rhosinAlpha - alpha * sinRhoCosAlpha)
            ) / (halfPi_RhoCosAlpha * halfPi_RhoCosAlpha));
        }
        r = rho + gingeryIntegrate(s_, rhoCosAlpha, x) * (pi$1 - rho) / gingeryIntegrate(s_, rhoCosAlpha, pi$1);
        theta = theta0 + alpha;
        x = r * cos$1(theta);
        y = r * sin$1(theta);
      }
      return azimuthalEquidistantRaw.invert(x, y);
    };

    return forward;
  }

  function gingeryLength(alpha, k) {
    return function(x) {
      var y_ = alpha * cos$1(x);
      if (x < halfPi$1) y_ -= k;
      return sqrt$1(1 + y_ * y_);
    };
  }

  // Numerical integration: trapezoidal rule.
  function gingeryIntegrate(f, a, b) {
    var n = 50,
        h = (b - a) / n,
        s = f(a) + f(b);
    for (var i = 1, x = a; i < n; ++i) s += 2 * f(x += h);
    return s * 0.5 * h;
  }

  function gingery() {
    var n = 6,
        rho = 30 * radians$1,
        cRho = cos$1(rho),
        sRho = sin$1(rho),
        m = projectionMutator(gingeryRaw),
        p = m(rho, n),
        stream_ = p.stream,
        epsilon = 1e-2,
        cr = -cos$1(epsilon * radians$1),
        sr = sin$1(epsilon * radians$1);

    p.radius = function(_) {
      if (!arguments.length) return rho * degrees$1;
      cRho = cos$1(rho = _ * radians$1);
      sRho = sin$1(rho);
      return m(rho, n);
    };

    p.lobes = function(_) {
      if (!arguments.length) return n;
      return m(rho, n = +_);
    };

    p.stream = function(stream) {
      var rotate = p.rotate(),
          rotateStream = stream_(stream),
          sphereStream = (p.rotate([0, 0]), stream_(stream));
      p.rotate(rotate);
      rotateStream.sphere = function() {
        sphereStream.polygonStart(), sphereStream.lineStart();
        for (var i = 0, delta = 2 * pi$1 / n, phi = 0; i < n; ++i, phi -= delta) {
          sphereStream.point(atan2$1(sr * cos$1(phi), cr) * degrees$1, asin$1(sr * sin$1(phi)) * degrees$1);
          sphereStream.point(atan2$1(sRho * cos$1(phi - delta / 2), cRho) * degrees$1, asin$1(sRho * sin$1(phi - delta / 2)) * degrees$1);
        }
        sphereStream.lineEnd(), sphereStream.polygonEnd();
      };
      return rotateStream;
    };

    return p
        .rotate([90, -40])
        .scale(91.7095)
        .clipAngle(180 - 1e-3);
  }

  function ginzburgPolyconicRaw(a, b, c, d, e, f, g, h) {
    if (arguments.length < 8) h = 0;

    function forward(lambda, phi) {
      if (!phi) return [a * lambda / pi$1, 0];
      var phi2 = phi * phi,
          xB = a + phi2 * (b + phi2 * (c + phi2 * d)),
          yB = phi * (e - 1 + phi2 * (f - h + phi2 * g)),
          m = (xB * xB + yB * yB) / (2 * yB),
          alpha = lambda * asin$1(xB / m) / pi$1;
      return [m * sin$1(alpha), phi * (1 + phi2 * h) + m * (1 - cos$1(alpha))];
    }

    forward.invert = function(x, y) {
      var lambda = pi$1 * x / a,
          phi = y,
          deltaLambda, deltaPhi, i = 50;
      do {
        var phi2 = phi * phi,
            xB = a + phi2 * (b + phi2 * (c + phi2 * d)),
            yB = phi * (e - 1 + phi2 * (f - h + phi2 * g)),
            p = xB * xB + yB * yB,
            q = 2 * yB,
            m = p / q,
            m2 = m * m,
            dAlphadLambda = asin$1(xB / m) / pi$1,
            alpha = lambda * dAlphadLambda,
            xB2 = xB * xB,
            dxBdPhi = (2 * b + phi2 * (4 * c + phi2 * 6 * d)) * phi,
            dyBdPhi = e + phi2 * (3 * f + phi2 * 5 * g),
            dpdPhi = 2 * (xB * dxBdPhi + yB * (dyBdPhi - 1)),
            dqdPhi = 2 * (dyBdPhi - 1),
            dmdPhi = (dpdPhi * q - p * dqdPhi) / (q * q),
            cosAlpha = cos$1(alpha),
            sinAlpha = sin$1(alpha),
            mcosAlpha = m * cosAlpha,
            msinAlpha = m * sinAlpha,
            dAlphadPhi = ((lambda / pi$1) * (1 / sqrt$1(1 - xB2 / m2)) * (dxBdPhi * m - xB * dmdPhi)) / m2,
            fx = msinAlpha - x,
            fy = phi * (1 + phi2 * h) + m - mcosAlpha - y,
            deltaxDeltaPhi = dmdPhi * sinAlpha + mcosAlpha * dAlphadPhi,
            deltaxDeltaLambda = mcosAlpha * dAlphadLambda,
            deltayDeltaPhi = 1 + dmdPhi - (dmdPhi * cosAlpha - msinAlpha * dAlphadPhi),
            deltayDeltaLambda = msinAlpha * dAlphadLambda,
            denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda;
        if (!denominator) break;
        lambda -= deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator;
        phi -= deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
      } while ((abs$1(deltaLambda) > epsilon$1 || abs$1(deltaPhi) > epsilon$1) && --i > 0);
      return [lambda, phi];
    };

    return forward;
  }

  var ginzburg4Raw = ginzburgPolyconicRaw(2.8284, -1.6988, 0.75432, -0.18071, 1.76003, -0.38914, 0.042555);

  function ginzburg4() {
    return projection(ginzburg4Raw)
        .scale(149.995);
  }

  var ginzburg5Raw = ginzburgPolyconicRaw(2.583819, -0.835827, 0.170354, -0.038094, 1.543313, -0.411435,0.082742);

  function ginzburg5() {
    return projection(ginzburg5Raw)
        .scale(153.93);
  }

  var ginzburg6Raw = ginzburgPolyconicRaw(5 / 6 * pi$1, -0.62636, -0.0344, 0, 1.3493, -0.05524, 0, 0.045);

  function ginzburg6() {
    return projection(ginzburg6Raw)
        .scale(130.945);
  }

  function ginzburg8Raw(lambda, phi) {
    var lambda2 = lambda * lambda,
        phi2 = phi * phi;
    return [
      lambda * (1 - 0.162388 * phi2) * (0.87 - 0.000952426 * lambda2 * lambda2),
      phi * (1 + phi2 / 12)
    ];
  }

  ginzburg8Raw.invert = function(x, y) {
    var lambda = x,
        phi = y,
        i = 50, delta;
    do {
      var phi2 = phi * phi;
      phi -= delta = (phi * (1 + phi2 / 12) - y) / (1 + phi2 / 4);
    } while (abs$1(delta) > epsilon$1 && --i > 0);
    i = 50;
    x /= 1 -0.162388 * phi2;
    do {
      var lambda4 = (lambda4 = lambda * lambda) * lambda4;
      lambda -= delta = (lambda * (0.87 - 0.000952426 * lambda4) - x) / (0.87 - 0.00476213 * lambda4);
    } while (abs$1(delta) > epsilon$1 && --i > 0);
    return [lambda, phi];
  };

  function ginzburg8() {
    return projection(ginzburg8Raw)
        .scale(131.747);
  }

  var ginzburg9Raw = ginzburgPolyconicRaw(2.6516, -0.76534, 0.19123, -0.047094, 1.36289, -0.13965,0.031762);

  function ginzburg9() {
    return projection(ginzburg9Raw)
        .scale(131.087);
  }

  function squareRaw(project) {
    var dx = project(halfPi$1, 0)[0] - project(-halfPi$1, 0)[0];

    function projectSquare(lambda, phi) {
      var s = lambda > 0 ? -0.5 : 0.5,
          point = project(lambda + s * pi$1, phi);
      point[0] -= s * dx;
      return point;
    }

    if (project.invert) projectSquare.invert = function(x, y) {
      var s = x > 0 ? -0.5 : 0.5,
          location = project.invert(x + s * dx, y),
          lambda = location[0] - s * pi$1;
      if (lambda < -pi$1) lambda += 2 * pi$1;
      else if (lambda > pi$1) lambda -= 2 * pi$1;
      location[0] = lambda;
      return location;
    };

    return projectSquare;
  }

  function gringortenRaw(lambda, phi) {
    var sLambda = sign$1(lambda),
        sPhi = sign$1(phi),
        cosPhi = cos$1(phi),
        x = cos$1(lambda) * cosPhi,
        y = sin$1(lambda) * cosPhi,
        z = sin$1(sPhi * phi);
    lambda = abs$1(atan2$1(y, z));
    phi = asin$1(x);
    if (abs$1(lambda - halfPi$1) > epsilon$1) lambda %= halfPi$1;
    var point = gringortenHexadecant(lambda > pi$1 / 4 ? halfPi$1 - lambda : lambda, phi);
    if (lambda > pi$1 / 4) z = point[0], point[0] = -point[1], point[1] = -z;
    return (point[0] *= sLambda, point[1] *= -sPhi, point);
  }

  gringortenRaw.invert = function(x, y) {
    if (abs$1(x) > 1) x = sign$1(x) * 2 - x;
    if (abs$1(y) > 1) y = sign$1(y) * 2 - y;
    var sx = sign$1(x),
        sy = sign$1(y),
        x0 = -sx * x,
        y0 = -sy * y,
        t = y0 / x0 < 1,
        p = gringortenHexadecantInvert(t ? y0 : x0, t ? x0 : y0),
        lambda = p[0],
        phi = p[1],
        cosPhi = cos$1(phi);
    if (t) lambda = -halfPi$1 - lambda;
    return [sx * (atan2$1(sin$1(lambda) * cosPhi, -sin$1(phi)) + pi$1), sy * asin$1(cos$1(lambda) * cosPhi)];
  };

  function gringortenHexadecant(lambda, phi) {
    if (phi === halfPi$1) return [0, 0];

    var sinPhi = sin$1(phi),
        r = sinPhi * sinPhi,
        r2 = r * r,
        j = 1 + r2,
        k = 1 + 3 * r2,
        q = 1 - r2,
        z = asin$1(1 / sqrt$1(j)),
        v = q + r * j * z,
        p2 = (1 - sinPhi) / v,
        p = sqrt$1(p2),
        a2 = p2 * j,
        a = sqrt$1(a2),
        h = p * q,
        x,
        i;

    if (lambda === 0) return [0, -(h + r * a)];

    var cosPhi = cos$1(phi),
        secPhi = 1 / cosPhi,
        drdPhi = 2 * sinPhi * cosPhi,
        dvdPhi = (-3 * r + z * k) * drdPhi,
        dp2dPhi = (-v * cosPhi - (1 - sinPhi) * dvdPhi) / (v * v),
        dpdPhi = (0.5 * dp2dPhi) / p,
        dhdPhi = q * dpdPhi - 2 * r * p * drdPhi,
        dra2dPhi = r * j * dp2dPhi + p2 * k * drdPhi,
        mu = -secPhi * drdPhi,
        nu = -secPhi * dra2dPhi,
        zeta = -2 * secPhi * dhdPhi,
        lambda1 = 4 * lambda / pi$1,
        delta;

    // Slower but accurate bisection method.
    if (lambda > 0.222 * pi$1 || phi < pi$1 / 4 && lambda > 0.175 * pi$1) {
      x = (h + r * sqrt$1(a2 * (1 + r2) - h * h)) / (1 + r2);
      if (lambda > pi$1 / 4) return [x, x];
      var x1 = x, x0 = 0.5 * x;
      x = 0.5 * (x0 + x1), i = 50;
      do {
        var g = sqrt$1(a2 - x * x),
            f = (x * (zeta + mu * g) + nu * asin$1(x / a)) - lambda1;
        if (!f) break;
        if (f < 0) x0 = x;
        else x1 = x;
        x = 0.5 * (x0 + x1);
      } while (abs$1(x1 - x0) > epsilon$1 && --i > 0);
    }

    // Newton-Raphson.
    else {
      x = epsilon$1, i = 25;
      do {
        var x2 = x * x,
            g2 = sqrt$1(a2 - x2),
            zetaMug = zeta + mu * g2,
            f2 = x * zetaMug + nu * asin$1(x / a) - lambda1,
            df = zetaMug + (nu - mu * x2) / g2;
        x -= delta = g2 ? f2 / df : 0;
      } while (abs$1(delta) > epsilon$1 && --i > 0);
    }

    return [x, -h - r * sqrt$1(a2 - x * x)];
  }

  function gringortenHexadecantInvert(x, y) {
    var x0 = 0,
        x1 = 1,
        r = 0.5,
        i = 50;

    while (true) {
      var r2 = r * r,
          sinPhi = sqrt$1(r),
          z = asin$1(1 / sqrt$1(1 + r2)),
          v = (1 - r2) + r * (1 + r2) * z,
          p2 = (1 - sinPhi) / v,
          p = sqrt$1(p2),
          a2 = p2 * (1 + r2),
          h = p * (1 - r2),
          g2 = a2 - x * x,
          g = sqrt$1(g2),
          y0 = y + h + r * g;
      if (abs$1(x1 - x0) < epsilon2$1 || --i === 0 || y0 === 0) break;
      if (y0 > 0) x0 = r;
      else x1 = r;
      r = 0.5 * (x0 + x1);
    }

    if (!i) return null;

    var phi = asin$1(sinPhi),
        cosPhi = cos$1(phi),
        secPhi = 1 / cosPhi,
        drdPhi = 2 * sinPhi * cosPhi,
        dvdPhi = (-3 * r + z * (1 + 3 * r2)) * drdPhi,
        dp2dPhi = (-v * cosPhi - (1 - sinPhi) * dvdPhi) / (v * v),
        dpdPhi = 0.5 * dp2dPhi / p,
        dhdPhi = (1 - r2) * dpdPhi - 2 * r * p * drdPhi,
        zeta = -2 * secPhi * dhdPhi,
        mu = -secPhi * drdPhi,
        nu = -secPhi * (r * (1 + r2) * dp2dPhi + p2 * (1 + 3 * r2) * drdPhi);

    return [pi$1 / 4 * (x * (zeta + mu * g) + nu * asin$1(x / sqrt$1(a2))), phi];
  }

  function gringorten() {
    return projection(squareRaw(gringortenRaw))
        .scale(239.75);
  }

  // Returns [sn, cn, dn](u + iv|m).
  function ellipticJi(u, v, m) {
    var a, b, c;
    if (!u) {
      b = ellipticJ(v, 1 - m);
      return [
        [0, b[0] / b[1]],
        [1 / b[1], 0],
        [b[2] / b[1], 0]
      ];
    }
    a = ellipticJ(u, m);
    if (!v) return [[a[0], 0], [a[1], 0], [a[2], 0]];
    b = ellipticJ(v, 1 - m);
    c = b[1] * b[1] + m * a[0] * a[0] * b[0] * b[0];
    return [
      [a[0] * b[2] / c, a[1] * a[2] * b[0] * b[1] / c],
      [a[1] * b[1] / c, -a[0] * a[2] * b[0] * b[2] / c],
      [a[2] * b[1] * b[2] / c, -m * a[0] * a[1] * b[0] / c]
    ];
  }

  // Returns [sn, cn, dn, ph](u|m).
  function ellipticJ(u, m) {
    var ai, b, phi, t, twon;
    if (m < epsilon$1) {
      t = sin$1(u);
      b = cos$1(u);
      ai = m * (u - t * b) / 4;
      return [
        t - ai * b,
        b + ai * t,
        1 - m * t * t / 2,
        u - ai
      ];
    }
    if (m >= 1 - epsilon$1) {
      ai = (1 - m) / 4;
      b = cosh(u);
      t = tanh(u);
      phi = 1 / b;
      twon = b * sinh(u);
      return [
        t + ai * (twon - u) / (b * b),
        phi - ai * t * phi * (twon - u),
        phi + ai * t * phi * (twon + u),
        2 * atan$1(exp$1(u)) - halfPi$1 + ai * (twon - u) / b
      ];
    }

    var a = [1, 0, 0, 0, 0, 0, 0, 0, 0],
        c = [sqrt$1(m), 0, 0, 0, 0, 0, 0, 0, 0],
        i = 0;
    b = sqrt$1(1 - m);
    twon = 1;

    while (abs$1(c[i] / a[i]) > epsilon$1 && i < 8) {
      ai = a[i++];
      c[i] = (ai - b) / 2;
      a[i] = (ai + b) / 2;
      b = sqrt$1(ai * b);
      twon *= 2;
    }

    phi = twon * a[i] * u;
    do {
      t = c[i] * sin$1(b = phi) / a[i];
      phi = (asin$1(t) + phi) / 2;
    } while (--i);

    return [sin$1(phi), t = cos$1(phi), t / cos$1(phi - b), phi];
  }

  // Calculate F(phi+iPsi|m).
  // See Abramowitz and Stegun, 17.4.11.
  function ellipticFi(phi, psi, m) {
    var r = abs$1(phi),
        i = abs$1(psi),
        sinhPsi = sinh(i);
    if (r) {
      var cscPhi = 1 / sin$1(r),
          cotPhi2 = 1 / (tan$1(r) * tan$1(r)),
          b = -(cotPhi2 + m * (sinhPsi * sinhPsi * cscPhi * cscPhi) - 1 + m),
          c = (m - 1) * cotPhi2,
          cotLambda2 = (-b + sqrt$1(b * b - 4 * c)) / 2;
      return [
        ellipticF(atan$1(1 / sqrt$1(cotLambda2)), m) * sign$1(phi),
        ellipticF(atan$1(sqrt$1((cotLambda2 / cotPhi2 - 1) / m)), 1 - m) * sign$1(psi)
      ];
    }
    return [
      0,
      ellipticF(atan$1(sinhPsi), 1 - m) * sign$1(psi)
    ];
  }

  // Calculate F(phi|m) where m = k² = sin²α.
  // See Abramowitz and Stegun, 17.6.7.
  function ellipticF(phi, m) {
    if (!m) return phi;
    if (m === 1) return log$1(tan$1(phi / 2 + quarterPi$1));
    var a = 1,
        b = sqrt$1(1 - m),
        c = sqrt$1(m);
    for (var i = 0; abs$1(c) > epsilon$1; i++) {
      if (phi % pi$1) {
        var dPhi = atan$1(b * tan$1(phi) / a);
        if (dPhi < 0) dPhi += pi$1;
        phi += dPhi + ~~(phi / pi$1) * pi$1;
      } else phi += phi;
      c = (a + b) / 2;
      b = sqrt$1(a * b);
      c = ((a = c) - b) / 2;
    }
    return phi / (pow$1(2, i) * a);
  }

  function guyouRaw(lambda, phi) {
    var k_ = (sqrt2 - 1) / (sqrt2 + 1),
        k = sqrt$1(1 - k_ * k_),
        K = ellipticF(halfPi$1, k * k),
        f = -1,
        psi = log$1(tan$1(pi$1 / 4 + abs$1(phi) / 2)),
        r = exp$1(f * psi) / sqrt$1(k_),
        at = guyouComplexAtan(r * cos$1(f * lambda), r * sin$1(f * lambda)),
        t = ellipticFi(at[0], at[1], k * k);
    return [-t[1], (phi >= 0 ? 1 : -1) * (0.5 * K - t[0])];
  }

  function guyouComplexAtan(x, y) {
    var x2 = x * x,
        y_1 = y + 1,
        t = 1 - x2 - y * y;
    return [
     0.5 * ((x >= 0 ? halfPi$1 : -halfPi$1) - atan2$1(t, 2 * x)),
      -0.25 * log$1(t * t + 4 * x2) +0.5 * log$1(y_1 * y_1 + x2)
    ];
  }

  function guyouComplexDivide(a, b) {
    var denominator = b[0] * b[0] + b[1] * b[1];
    return [
      (a[0] * b[0] + a[1] * b[1]) / denominator,
      (a[1] * b[0] - a[0] * b[1]) / denominator
    ];
  }

  guyouRaw.invert = function(x, y) {
    var k_ = (sqrt2 - 1) / (sqrt2 + 1),
        k = sqrt$1(1 - k_ * k_),
        K = ellipticF(halfPi$1, k * k),
        f = -1,
        j = ellipticJi(0.5 * K - y, -x, k * k),
        tn = guyouComplexDivide(j[0], j[1]),
        lambda = atan2$1(tn[1], tn[0]) / f;
    return [
      lambda,
      2 * atan$1(exp$1(0.5 / f * log$1(k_ * tn[0] * tn[0] + k_ * tn[1] * tn[1]))) - halfPi$1
    ];
  };

  function guyou() {
    return projection(squareRaw(guyouRaw))
        .scale(151.496);
  }

  function hammerRetroazimuthalRaw(phi0) {
    var sinPhi0 = sin$1(phi0),
        cosPhi0 = cos$1(phi0),
        rotate = hammerRetroazimuthalRotation(phi0);

    rotate.invert = hammerRetroazimuthalRotation(-phi0);

    function forward(lambda, phi) {
      var p = rotate(lambda, phi);
      lambda = p[0], phi = p[1];
      var sinPhi = sin$1(phi),
          cosPhi = cos$1(phi),
          cosLambda = cos$1(lambda),
          z = acos$1(sinPhi0 * sinPhi + cosPhi0 * cosPhi * cosLambda),
          sinz = sin$1(z),
          K = abs$1(sinz) > epsilon$1 ? z / sinz : 1;
      return [
        K * cosPhi0 * sin$1(lambda),
        (abs$1(lambda) > halfPi$1 ? K : -K) // rotate for back hemisphere
          * (sinPhi0 * cosPhi - cosPhi0 * sinPhi * cosLambda)
      ];
    }

    forward.invert = function(x, y) {
      var rho = sqrt$1(x * x + y * y),
          sinz = -sin$1(rho),
          cosz = cos$1(rho),
          a = rho * cosz,
          b = -y * sinz,
          c = rho * sinPhi0,
          d = sqrt$1(a * a + b * b - c * c),
          phi = atan2$1(a * c + b * d, b * c - a * d),
          lambda = (rho > halfPi$1 ? -1 : 1) * atan2$1(x * sinz, rho * cos$1(phi) * cosz + y * sin$1(phi) * sinz);
      return rotate.invert(lambda, phi);
    };

    return forward;
  }

  // Latitudinal rotation by phi0.
  // Temporary hack until D3 supports arbitrary small-circle clipping origins.
  function hammerRetroazimuthalRotation(phi0) {
    var sinPhi0 = sin$1(phi0),
        cosPhi0 = cos$1(phi0);

    return function(lambda, phi) {
      var cosPhi = cos$1(phi),
          x = cos$1(lambda) * cosPhi,
          y = sin$1(lambda) * cosPhi,
          z = sin$1(phi);
      return [
        atan2$1(y, x * cosPhi0 - z * sinPhi0),
        asin$1(z * cosPhi0 + x * sinPhi0)
      ];
    };
  }

  function hammerRetroazimuthal() {
    var phi0 = 0,
        m = projectionMutator(hammerRetroazimuthalRaw),
        p = m(phi0),
        rotate_ = p.rotate,
        stream_ = p.stream,
        circle = geoCircle();

    p.parallel = function(_) {
      if (!arguments.length) return phi0 * degrees$1;
      var r = p.rotate();
      return m(phi0 = _ * radians$1).rotate(r);
    };

    // Temporary hack; see hammerRetroazimuthalRotation.
    p.rotate = function(_) {
      if (!arguments.length) return (_ = rotate_.call(p), _[1] += phi0 * degrees$1, _);
      rotate_.call(p, [_[0], _[1] - phi0 * degrees$1]);
      circle.center([-_[0], -_[1]]);
      return p;
    };

    p.stream = function(stream) {
      stream = stream_(stream);
      stream.sphere = function() {
        stream.polygonStart();
        var epsilon = 1e-2,
            ring = circle.radius(90 - epsilon)().coordinates[0],
            n = ring.length - 1,
            i = -1,
            p;
        stream.lineStart();
        while (++i < n) stream.point((p = ring[i])[0], p[1]);
        stream.lineEnd();
        ring = circle.radius(90 + epsilon)().coordinates[0];
        n = ring.length - 1;
        stream.lineStart();
        while (--i >= 0) stream.point((p = ring[i])[0], p[1]);
        stream.lineEnd();
        stream.polygonEnd();
      };
      return stream;
    };

    return p
        .scale(79.4187)
        .parallel(45)
        .clipAngle(180 - 1e-3);
  }

  var healpixParallel = 41 + 48 / 36 + 37 / 3600, // for K=3; TODO automate
      healpixLambert = cylindricalEqualAreaRaw$1(0);

  function healpixRaw(H) {
    var phi0 = healpixParallel * radians$1,
        dx = collignonRaw(pi$1, phi0)[0] - collignonRaw(-pi$1, phi0)[0],
        y0 = healpixLambert(0, phi0)[1],
        y1 = collignonRaw(0, phi0)[1],
        dy1 = sqrtPi - y1,
        k = tau$1 / H,
        w = 4 / tau$1,
        h = y0 + (dy1 * dy1 * 4) / tau$1;

    function forward(lambda, phi) {
      var point,
          phi2 = abs$1(phi);
      if (phi2 > phi0) {
        var i = min$1(H - 1, max$1(0, floor$1((lambda + pi$1) / k)));
        lambda += pi$1 * (H - 1) / H - i * k;
        point = collignonRaw(lambda, phi2);
        point[0] = point[0] * tau$1 / dx - tau$1 * (H - 1) / (2 * H) + i * tau$1 / H;
        point[1] = y0 + (point[1] - y1) * 4 * dy1 / tau$1;
        if (phi < 0) point[1] = -point[1];
      } else {
        point = healpixLambert(lambda, phi);
      }
      point[0] *= w, point[1] /= h;
      return point;
    }

    forward.invert = function(x, y) {
      x /= w, y *= h;
      var y2 = abs$1(y);
      if (y2 > y0) {
        var i = min$1(H - 1, max$1(0, floor$1((x + pi$1) / k)));
        x = (x + pi$1 * (H - 1) / H - i * k) * dx / tau$1;
        var point = collignonRaw.invert(x, 0.25 * (y2 - y0) * tau$1 / dy1 + y1);
        point[0] -= pi$1 * (H - 1) / H - i * k;
        if (y < 0) point[1] = -point[1];
        return point;
      }
      return healpixLambert.invert(x, y);
    };

    return forward;
  }

  function sphere(step) {
    return {
      type: "Polygon",
      coordinates: [
        range$1(-180, 180 + step / 2, step).map(function(x, i) { return [x, i & 1 ? 90 - 1e-6 : healpixParallel]; })
        .concat(range$1(180, -180 - step / 2, -step).map(function(x, i) { return [x, i & 1 ? -90 + 1e-6 : -healpixParallel]; }))
      ]
    };
  }

  function healpix() {
    var H = 4,
        m = projectionMutator(healpixRaw),
        p = m(H),
        stream_ = p.stream;

    p.lobes = function(_) {
      return arguments.length ? m(H = +_) : H;
    };

    p.stream = function(stream) {
      var rotate = p.rotate(),
          rotateStream = stream_(stream),
          sphereStream = (p.rotate([0, 0]), stream_(stream));
      p.rotate(rotate);
      rotateStream.sphere = function() { geoStream(sphere(180 / H), sphereStream); };
      return rotateStream;
    };

    return p
        .scale(239.75);
  }

  function hillRaw(K) {
    var L = 1 + K,
        sinBt = sin$1(1 / L),
        Bt = asin$1(sinBt),
        A = 2 * sqrt$1(pi$1 / (B = pi$1 + 4 * Bt * L)),
        B,
        rho0 = 0.5 * A * (L + sqrt$1(K * (2 + K))),
        K2 = K * K,
        L2 = L * L;

    function forward(lambda, phi) {
      var t = 1 - sin$1(phi),
          rho,
          omega;
      if (t && t < 2) {
        var theta = halfPi$1 - phi, i = 25, delta;
        do {
          var sinTheta = sin$1(theta),
              cosTheta = cos$1(theta),
              Bt_Bt1 = Bt + atan2$1(sinTheta, L - cosTheta),
              C = 1 + L2 - 2 * L * cosTheta;
          theta -= delta = (theta - K2 * Bt - L * sinTheta + C * Bt_Bt1 -0.5 * t * B) / (2 * L * sinTheta * Bt_Bt1);
        } while (abs$1(delta) > epsilon2$1 && --i > 0);
        rho = A * sqrt$1(C);
        omega = lambda * Bt_Bt1 / pi$1;
      } else {
        rho = A * (K + t);
        omega = lambda * Bt / pi$1;
      }
      return [
        rho * sin$1(omega),
        rho0 - rho * cos$1(omega)
      ];
    }

    forward.invert = function(x, y) {
      var rho2 = x * x + (y -= rho0) * y,
          cosTheta = (1 + L2 - rho2 / (A * A)) / (2 * L),
          theta = acos$1(cosTheta),
          sinTheta = sin$1(theta),
          Bt_Bt1 = Bt + atan2$1(sinTheta, L - cosTheta);
      return [
        asin$1(x / sqrt$1(rho2)) * pi$1 / Bt_Bt1,
        asin$1(1 - 2 * (theta - K2 * Bt - L * sinTheta + (1 + L2 - 2 * L * cosTheta) * Bt_Bt1) / B)
      ];
    };

    return forward;
  }

  function hill() {
    var K = 1,
        m = projectionMutator(hillRaw),
        p = m(K);

    p.ratio = function(_) {
      return arguments.length ? m(K = +_) : K;
    };

    return p
        .scale(167.774)
        .center([0, 18.67]);
  }

  var sinuMollweidePhi = 0.7109889596207567;

  var sinuMollweideY = 0.0528035274542;

  function sinuMollweideRaw(lambda, phi) {
    return phi > -sinuMollweidePhi
        ? (lambda = mollweideRaw(lambda, phi), lambda[1] += sinuMollweideY, lambda)
        : sinusoidalRaw(lambda, phi);
  }

  sinuMollweideRaw.invert = function(x, y) {
    return y > -sinuMollweidePhi
        ? mollweideRaw.invert(x, y - sinuMollweideY)
        : sinusoidalRaw.invert(x, y);
  };

  function sinuMollweide() {
    return projection(sinuMollweideRaw)
        .rotate([-20, -55])
        .scale(164.263)
        .center([0, -5.4036]);
  }

  function homolosineRaw(lambda, phi) {
    return abs$1(phi) > sinuMollweidePhi
        ? (lambda = mollweideRaw(lambda, phi), lambda[1] -= phi > 0 ? sinuMollweideY : -sinuMollweideY, lambda)
        : sinusoidalRaw(lambda, phi);
  }

  homolosineRaw.invert = function(x, y) {
    return abs$1(y) > sinuMollweidePhi
        ? mollweideRaw.invert(x, y + (y > 0 ? sinuMollweideY : -sinuMollweideY))
        : sinusoidalRaw.invert(x, y);
  };

  function homolosine() {
    return projection(homolosineRaw)
        .scale(152.63);
  }

  // https://github.com/scijs/integrate-adaptive-simpson

  // This algorithm adapted from pseudocode in:
  // http://www.math.utk.edu/~ccollins/refs/Handouts/rich.pdf
  function adsimp (f, a, b, fa, fm, fb, V0, tol, maxdepth, depth, state) {
    if (state.nanEncountered) {
      return NaN;
    }

    var h, f1, f2, sl, sr, s2, m, V1, V2, err;

    h = b - a;
    f1 = f(a + h * 0.25);
    f2 = f(b - h * 0.25);

    // Simple check for NaN:
    if (isNaN(f1)) {
      state.nanEncountered = true;
      return;
    }

    // Simple check for NaN:
    if (isNaN(f2)) {
      state.nanEncountered = true;
      return;
    }

    sl = h * (fa + 4 * f1 + fm) / 12;
    sr = h * (fm + 4 * f2 + fb) / 12;
    s2 = sl + sr;
    err = (s2 - V0) / 15;

    if (depth > maxdepth) {
      state.maxDepthCount++;
      return s2 + err;
    } else if (Math.abs(err) < tol) {
      return s2 + err;
    } else {
      m = a + h * 0.5;

      V1 = adsimp(f, a, m, fa, f1, fm, sl, tol * 0.5, maxdepth, depth + 1, state);

      if (isNaN(V1)) {
        state.nanEncountered = true;
        return NaN;
      }

      V2 = adsimp(f, m, b, fm, f2, fb, sr, tol * 0.5, maxdepth, depth + 1, state);

      if (isNaN(V2)) {
        state.nanEncountered = true;
        return NaN;
      }

      return V1 + V2;
    }
  }

  function integrate (f, a, b, tol, maxdepth) {
    var state = {
      maxDepthCount: 0,
      nanEncountered: false
    };

    if (tol === undefined) {
      tol = 1e-8;
    }
    if (maxdepth === undefined) {
      maxdepth = 20;
    }

    var fa = f(a);
    var fm = f(0.5 * (a + b));
    var fb = f(b);

    var V0 = (fa + 4 * fm + fb) * (b - a) / 6;

    var result = adsimp(f, a, b, fa, fm, fb, V0, tol, maxdepth, 1, state);

  /*
    if (state.maxDepthCount > 0 && console && console.warn) {
      console.warn('integrate-adaptive-simpson: Warning: maximum recursion depth (' + maxdepth + ') reached ' + state.maxDepthCount + ' times');
    }

    if (state.nanEncountered && console && console.warn) {
      console.warn('integrate-adaptive-simpson: Warning: NaN encountered. Halting early.');
    }
  */

    return result;
  }

  function hyperellipticalRaw(alpha, k, gamma) {

    function elliptic (f) {
      return alpha + (1 - alpha) * pow$1(1 - pow$1(f, k), 1 / k);
    }

    function z(f) {
      return integrate(elliptic, 0, f, 1e-4);
    }

    var G = 1 / z(1),
        n = 1000,
        m = (1 + 1e-8) * G,
        approx = [];
    for (var i = 0; i <= n; i++)
        approx.push(z(i / n) * m);

    function Y(sinphi) {
      var rmin = 0, rmax = n, r = n >> 1;
      do {
        if (approx[r] > sinphi) rmax = r; else rmin = r;
        r = (rmin + rmax) >> 1;
      } while (r > rmin);
      var u = approx[r + 1] - approx[r];
      if (u) u = (sinphi - approx[r + 1]) / u;
      return (r + 1 + u) / n;
    }

    var ratio = 2 * Y(1) / pi$1 * G / gamma;

    var forward = function(lambda, phi) {
      var y = Y(abs$1(sin$1(phi))),
          x = elliptic(y) * lambda;
      y /= ratio;
      return [ x, (phi >= 0) ? y : -y ];
    };

    forward.invert = function(x, y) {
      var phi;
      y *= ratio;
      if (abs$1(y) < 1) phi = sign$1(y) * asin$1(z(abs$1(y)) * G);
      return [ x / elliptic(abs$1(y)), phi ];
    };

    return forward;
  }

  function hyperelliptical() {
    var alpha = 0,
        k = 2.5,
        gamma = 1.183136, // affine = sqrt(2 * gamma / pi) = 0.8679
        m = projectionMutator(hyperellipticalRaw),
        p = m(alpha, k, gamma);

    p.alpha = function(_) {
      return arguments.length ? m(alpha = +_, k, gamma) : alpha;
    };

    p.k = function(_) {
      return arguments.length ? m(alpha, k = +_, gamma) : k;
    };

    p.gamma = function(_) {
      return arguments.length ? m(alpha, k, gamma = +_) : gamma;
    };

    return p
        .scale(152.63);
  }

  function pointEqual$1(a, b) {
    return abs$1(a[0] - b[0]) < epsilon$1 && abs$1(a[1] - b[1]) < epsilon$1;
  }

  function interpolateLine(coordinates, m) {
    var i = -1,
        n = coordinates.length,
        p0 = coordinates[0],
        p1,
        dx,
        dy,
        resampled = [];
    while (++i < n) {
      p1 = coordinates[i];
      dx = (p1[0] - p0[0]) / m;
      dy = (p1[1] - p0[1]) / m;
      for (var j = 0; j < m; ++j) resampled.push([p0[0] + j * dx, p0[1] + j * dy]);
      p0 = p1;
    }
    resampled.push(p1);
    return resampled;
  }

  function interpolateSphere(lobes) {
    var coordinates = [],
        lobe,
        lambda0, phi0, phi1,
        lambda2, phi2,
        i, n = lobes[0].length;

    // Northern Hemisphere
    for (i = 0; i < n; ++i) {
      lobe = lobes[0][i];
      lambda0 = lobe[0][0], phi0 = lobe[0][1], phi1 = lobe[1][1];
      lambda2 = lobe[2][0], phi2 = lobe[2][1];
      coordinates.push(interpolateLine([
        [lambda0 + epsilon$1, phi0 + epsilon$1],
        [lambda0 + epsilon$1, phi1 - epsilon$1],
        [lambda2 - epsilon$1, phi1 - epsilon$1],
        [lambda2 - epsilon$1, phi2 + epsilon$1]
      ], 30));
    }

    // Southern Hemisphere
    for (i = lobes[1].length - 1; i >= 0; --i) {
      lobe = lobes[1][i];
      lambda0 = lobe[0][0], phi0 = lobe[0][1], phi1 = lobe[1][1];
      lambda2 = lobe[2][0], phi2 = lobe[2][1];
      coordinates.push(interpolateLine([
        [lambda2 - epsilon$1, phi2 - epsilon$1],
        [lambda2 - epsilon$1, phi1 + epsilon$1],
        [lambda0 + epsilon$1, phi1 + epsilon$1],
        [lambda0 + epsilon$1, phi0 - epsilon$1]
      ], 30));
    }

    return {
      type: "Polygon",
      coordinates: [merge(coordinates)]
    };
  }

  function interrupt(project, lobes) {
    var sphere, bounds$$1;

    function forward(lambda, phi) {
      var sign = phi < 0 ? -1 : +1, lobe = lobes[+(phi < 0)];
      for (var i = 0, n = lobe.length - 1; i < n && lambda > lobe[i][2][0]; ++i);
      var p = project(lambda - lobe[i][1][0], phi);
      p[0] += project(lobe[i][1][0], sign * phi > sign * lobe[i][0][1] ? lobe[i][0][1] : phi)[0];
      return p;
    }

    // Assumes mutually exclusive bounding boxes for lobes.
    if (project.invert) forward.invert = function(x, y) {
      var bound = bounds$$1[+(y < 0)], lobe = lobes[+(y < 0)];
      for (var i = 0, n = bound.length; i < n; ++i) {
        var b = bound[i];
        if (b[0][0] <= x && x < b[1][0] && b[0][1] <= y && y < b[1][1]) {
          var p = project.invert(x - project(lobe[i][1][0], 0)[0], y);
          p[0] += lobe[i][1][0];
          return pointEqual$1(forward(p[0], p[1]), [x, y]) ? p : null;
        }
      }
    };

    var p = projection(forward),
        stream_ = p.stream;

    p.stream = function(stream) {
      var rotate = p.rotate(),
          rotateStream = stream_(stream),
          sphereStream = (p.rotate([0, 0]), stream_(stream));
      p.rotate(rotate);
      rotateStream.sphere = function() { geoStream(sphere, sphereStream); };
      return rotateStream;
    };
    
    p.lobes = function(_) {
      if (!arguments.length) return lobes.map(function(lobe) {
        return lobe.map(function(l) {
          return [
            [l[0][0] * degrees$1, l[0][1] * degrees$1],
            [l[1][0] * degrees$1, l[1][1] * degrees$1],
            [l[2][0] * degrees$1, l[2][1] * degrees$1]
          ];
        });
      });

      sphere = interpolateSphere(_);

      lobes = _.map(function(lobe) {
        return lobe.map(function(l) {
          return [
            [l[0][0] * radians$1, l[0][1] * radians$1],
            [l[1][0] * radians$1, l[1][1] * radians$1],
            [l[2][0] * radians$1, l[2][1] * radians$1]
          ];
        });
      });

      bounds$$1 = lobes.map(function(lobe) {
        return lobe.map(function(l) {
          var x0 = project(l[0][0], l[0][1])[0],
              x1 = project(l[2][0], l[2][1])[0],
              y0 = project(l[1][0], l[0][1])[1],
              y1 = project(l[1][0], l[1][1])[1],
              t;
          if (y0 > y1) t = y0, y0 = y1, y1 = t;
          return [[x0, y0], [x1, y1]];
        });
      });

      return p;
    };

    if (lobes != null) p.lobes(lobes);

    return p;
  }

  var lobes = [[ // northern hemisphere
    [[-180,   0], [-100,  90], [ -40,   0]],
    [[ -40,   0], [  30,  90], [ 180,   0]]
  ], [ // southern hemisphere
    [[-180,   0], [-160, -90], [-100,   0]],
    [[-100,   0], [ -60, -90], [ -20,   0]],
    [[ -20,   0], [  20, -90], [  80,   0]],
    [[  80,   0], [ 140, -90], [ 180,   0]]
  ]];

  function boggs$1() {
    return interrupt(boggsRaw, lobes)
        .scale(160.857);
  }

  var lobes$1 = [[ // northern hemisphere
    [[-180,   0], [-100,  90], [ -40,   0]],
    [[ -40,   0], [  30,  90], [ 180,   0]]
  ], [ // southern hemisphere
    [[-180,   0], [-160, -90], [-100,   0]],
    [[-100,   0], [ -60, -90], [ -20,   0]],
    [[ -20,   0], [  20, -90], [  80,   0]],
    [[  80,   0], [ 140, -90], [ 180,   0]]
  ]];

  function homolosine$1() {
    return interrupt(homolosineRaw, lobes$1)
        .scale(152.63);
  }

  var lobes$2 = [[ // northern hemisphere
    [[-180,   0], [-100,  90], [ -40,   0]],
    [[ -40,   0], [  30,  90], [ 180,   0]]
  ], [ // southern hemisphere
    [[-180,   0], [-160, -90], [-100,   0]],
    [[-100,   0], [ -60, -90], [ -20,   0]],
    [[ -20,   0], [  20, -90], [  80,   0]],
    [[  80,   0], [ 140, -90], [ 180,   0]]
  ]];

  function mollweide$1() {
    return interrupt(mollweideRaw, lobes$2)
        .scale(169.529);
  }

  var lobes$3 = [[ // northern hemisphere
    [[-180,   0], [ -90,  90], [   0,   0]],
    [[   0,   0], [  90,  90], [ 180,   0]]
  ], [ // southern hemisphere
    [[-180,   0], [ -90, -90], [   0,   0]],
    [[   0,   0], [  90, -90], [ 180,   0]]
  ]];

  function mollweideHemispheres() {
    return interrupt(mollweideRaw, lobes$3)
        .scale(169.529)
        .rotate([20, 0]);
  }

  var lobes$4 = [[ // northern hemisphere
    [[-180,  35], [ -30,  90], [   0,  35]],
    [[   0,  35], [  30,  90], [ 180,  35]]
  ], [ // southern hemisphere
    [[-180, -10], [-102, -90], [ -65, -10]],
    [[ -65, -10], [   5, -90], [  77, -10]],
    [[  77, -10], [ 103, -90], [ 180, -10]]
  ]];

  function sinuMollweide$1() {
    return interrupt(sinuMollweideRaw, lobes$4)
        .rotate([-20, -55])
        .scale(164.263)
        .center([0, -5.4036]);
  }

  var lobes$5 = [[ // northern hemisphere
    [[-180,   0], [-110,  90], [ -40,   0]],
    [[ -40,   0], [   0,  90], [  40,   0]],
    [[  40,   0], [ 110,  90], [ 180,   0]]
  ], [ // southern hemisphere
    [[-180,   0], [-110, -90], [ -40,   0]],
    [[ -40,   0], [   0, -90], [  40,   0]],
    [[  40,   0], [ 110, -90], [ 180,   0]]
  ]];

  function sinusoidal$1() {
    return interrupt(sinusoidalRaw, lobes$5)
        .scale(152.63)
        .rotate([-20, 0]);
  }

  function kavrayskiy7Raw(lambda, phi) {
    return [3 / tau$1 * lambda * sqrt$1(pi$1 * pi$1 / 3 - phi * phi), phi];
  }

  kavrayskiy7Raw.invert = function(x, y) {
    return [tau$1 / 3 * x / sqrt$1(pi$1 * pi$1 / 3 - y * y), y];
  };

  function kavrayskiy7() {
    return projection(kavrayskiy7Raw)
        .scale(158.837);
  }

  function lagrangeRaw(n) {

    function forward(lambda, phi) {
      if (abs$1(abs$1(phi) - halfPi$1) < epsilon$1) return [0, phi < 0 ? -2 : 2];
      var sinPhi = sin$1(phi),
          v = pow$1((1 + sinPhi) / (1 - sinPhi), n / 2),
          c = 0.5 * (v + 1 / v) + cos$1(lambda *= n);
      return [
        2 * sin$1(lambda) / c,
        (v - 1 / v) / c
      ];
    }

    forward.invert = function(x, y) {
      var y0 = abs$1(y);
      if (abs$1(y0 - 2) < epsilon$1) return x ? null : [0, sign$1(y) * halfPi$1];
      if (y0 > 2) return null;

      x /= 2, y /= 2;
      var x2 = x * x,
          y2 = y * y,
          t = 2 * y / (1 + x2 + y2); // tanh(nPhi)
      t = pow$1((1 + t) / (1 - t), 1 / n);
      return [
        atan2$1(2 * x, 1 - x2 - y2) / n,
        asin$1((t - 1) / (t + 1))
      ];
    };

    return forward;
  }

  function lagrange() {
    var n = 0.5,
        m = projectionMutator(lagrangeRaw),
        p = m(n);

    p.spacing = function(_) {
      return arguments.length ? m(n = +_) : n;
    };

    return p
        .scale(124.75);
  }

  var pi_sqrt2 = pi$1 / sqrt2;

  function larriveeRaw(lambda, phi) {
    return [
      lambda * (1 + sqrt$1(cos$1(phi))) / 2,
      phi / (cos$1(phi / 2) * cos$1(lambda / 6))
    ];
  }

  larriveeRaw.invert = function(x, y) {
    var x0 = abs$1(x),
        y0 = abs$1(y),
        lambda = epsilon$1,
        phi = halfPi$1;
    if (y0 < pi_sqrt2) phi *= y0 / pi_sqrt2;
    else lambda += 6 * acos$1(pi_sqrt2 / y0);
    for (var i = 0; i < 25; i++) {
      var sinPhi = sin$1(phi),
          sqrtcosPhi = sqrt$1(cos$1(phi)),
          sinPhi_2 = sin$1(phi / 2),
          cosPhi_2 = cos$1(phi / 2),
          sinLambda_6 = sin$1(lambda / 6),
          cosLambda_6 = cos$1(lambda / 6),
          f0 = 0.5 * lambda * (1 + sqrtcosPhi) - x0,
          f1 = phi / (cosPhi_2 * cosLambda_6) - y0,
          df0dPhi = sqrtcosPhi ? -0.25 * lambda * sinPhi / sqrtcosPhi : 0,
          df0dLambda = 0.5 * (1 + sqrtcosPhi),
          df1dPhi = (1 +0.5 * phi * sinPhi_2 / cosPhi_2) / (cosPhi_2 * cosLambda_6),
          df1dLambda = (phi / cosPhi_2) * (sinLambda_6 / 6) / (cosLambda_6 * cosLambda_6),
          denom = df0dPhi * df1dLambda - df1dPhi * df0dLambda,
          dPhi = (f0 * df1dLambda - f1 * df0dLambda) / denom,
          dLambda = (f1 * df0dPhi - f0 * df1dPhi) / denom;
      phi -= dPhi;
      lambda -= dLambda;
      if (abs$1(dPhi) < epsilon$1 && abs$1(dLambda) < epsilon$1) break;
    }
    return [x < 0 ? -lambda : lambda, y < 0 ? -phi : phi];
  };

  function larrivee() {
    return projection(larriveeRaw)
        .scale(97.2672);
  }

  function laskowskiRaw(lambda, phi) {
    var lambda2 = lambda * lambda, phi2 = phi * phi;
    return [
      lambda * (0.975534 + phi2 * (-0.119161 + lambda2 * -0.0143059 + phi2 * -0.0547009)),
      phi * (1.00384 + lambda2 * (0.0802894 + phi2 * -0.02855 + lambda2 * 0.000199025) + phi2 * (0.0998909 + phi2 * -0.0491032))
    ];
  }

  laskowskiRaw.invert = function(x, y) {
    var lambda = sign$1(x) * pi$1,
        phi = y / 2,
        i = 50;
    do {
      var lambda2 = lambda * lambda,
          phi2 = phi * phi,
          lambdaPhi = lambda * phi,
          fx = lambda * (0.975534 + phi2 * (-0.119161 + lambda2 * -0.0143059 + phi2 * -0.0547009)) - x,
          fy = phi * (1.00384 + lambda2 * (0.0802894 + phi2 * -0.02855 + lambda2 * 0.000199025) + phi2 * (0.0998909 + phi2 * -0.0491032)) - y,
          deltaxDeltaLambda = 0.975534 - phi2 * (0.119161 + 3 * lambda2 * 0.0143059 + phi2 * 0.0547009),
          deltaxDeltaPhi = -lambdaPhi * (2 * 0.119161 + 4 * 0.0547009 * phi2 + 2 * 0.0143059 * lambda2),
          deltayDeltaLambda = lambdaPhi * (2 * 0.0802894 + 4 * 0.000199025 * lambda2 + 2 * -0.02855 * phi2),
          deltayDeltaPhi = 1.00384 + lambda2 * (0.0802894 + 0.000199025 * lambda2) + phi2 * (3 * (0.0998909 - 0.02855 * lambda2) - 5 * 0.0491032 * phi2),
          denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda,
          deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator,
          deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
      lambda -= deltaLambda, phi -= deltaPhi;
    } while ((abs$1(deltaLambda) > epsilon$1 || abs$1(deltaPhi) > epsilon$1) && --i > 0);
    return i && [lambda, phi];
  };

  function laskowski() {
    return projection(laskowskiRaw)
        .scale(139.98);
  }

  function littrowRaw(lambda, phi) {
    return [
      sin$1(lambda) / cos$1(phi),
      tan$1(phi) * cos$1(lambda)
    ];
  }

  littrowRaw.invert = function(x, y) {
    var x2 = x * x,
        y2 = y * y,
        y2_1 = y2 + 1,
        x2_y2_1 = x2 + y2_1,
        cosPhi = x
            ? sqrt1_2 * sqrt$1((x2_y2_1 - sqrt$1(x2_y2_1 * x2_y2_1 - 4 * x2)) / x2)
            : 1 / sqrt$1(y2_1);
    return [
      asin$1(x * cosPhi),
      sign$1(y) * acos$1(cosPhi)
    ];
  };

  function littrow() {
    return projection(littrowRaw)
        .scale(144.049)
        .clipAngle(90 - 1e-3);
  }

  function loximuthalRaw(phi0) {
    var cosPhi0 = cos$1(phi0),
        tanPhi0 = tan$1(quarterPi$1 + phi0 / 2);

    function forward(lambda, phi) {
      var y = phi - phi0,
          x = abs$1(y) < epsilon$1 ? lambda * cosPhi0
              : abs$1(x = quarterPi$1 + phi / 2) < epsilon$1 || abs$1(abs$1(x) - halfPi$1) < epsilon$1
              ? 0 : lambda * y / log$1(tan$1(x) / tanPhi0);
      return [x, y];
    }

    forward.invert = function(x, y) {
      var lambda,
          phi = y + phi0;
      return [
        abs$1(y) < epsilon$1 ? x / cosPhi0
            : (abs$1(lambda = quarterPi$1 + phi / 2) < epsilon$1 || abs$1(abs$1(lambda) - halfPi$1) < epsilon$1) ? 0
            : x * log$1(tan$1(lambda) / tanPhi0) / y,
        phi
      ];
    };

    return forward;
  }

  function loximuthal() {
    return parallel1(loximuthalRaw)
        .parallel(40)
        .scale(158.837);
  }

  function millerRaw(lambda, phi) {
    return [lambda, 1.25 * log$1(tan$1(quarterPi$1 + 0.4 * phi))];
  }

  millerRaw.invert = function(x, y) {
    return [x, 2.5 * atan$1(exp$1(0.8 * y)) - 0.625 * pi$1];
  };

  function miller() {
    return projection(millerRaw)
        .scale(108.318);
  }

  function modifiedStereographicRaw(C) {
    var m = C.length - 1;

    function forward(lambda, phi) {
      var cosPhi = cos$1(phi),
          k = 2 / (1 + cosPhi * cos$1(lambda)),
          zr = k * cosPhi * sin$1(lambda),
          zi = k * sin$1(phi),
          i = m,
          w = C[i],
          ar = w[0],
          ai = w[1],
          t;
      while (--i >= 0) {
        w = C[i];
        ar = w[0] + zr * (t = ar) - zi * ai;
        ai = w[1] + zr * ai + zi * t;
      }
      ar = zr * (t = ar) - zi * ai;
      ai = zr * ai + zi * t;
      return [ar, ai];
    }

    forward.invert = function(x, y) {
      var i = 20,
          zr = x,
          zi = y;
      do {
        var j = m,
            w = C[j],
            ar = w[0],
            ai = w[1],
            br = 0,
            bi = 0,
            t;

        while (--j >= 0) {
          w = C[j];
          br = ar + zr * (t = br) - zi * bi;
          bi = ai + zr * bi + zi * t;
          ar = w[0] + zr * (t = ar) - zi * ai;
          ai = w[1] + zr * ai + zi * t;
        }
        br = ar + zr * (t = br) - zi * bi;
        bi = ai + zr * bi + zi * t;
        ar = zr * (t = ar) - zi * ai - x;
        ai = zr * ai + zi * t - y;

        var denominator = br * br + bi * bi, deltar, deltai;
        zr -= deltar = (ar * br + ai * bi) / denominator;
        zi -= deltai = (ai * br - ar * bi) / denominator;
      } while (abs$1(deltar) + abs$1(deltai) > epsilon$1 * epsilon$1 && --i > 0);

      if (i) {
        var rho = sqrt$1(zr * zr + zi * zi),
            c = 2 * atan$1(rho * 0.5),
            sinc = sin$1(c);
        return [atan2$1(zr * sinc, rho * cos$1(c)), rho ? asin$1(zi * sinc / rho) : 0];
      }
    };

    return forward;
  }

  var alaska = [[0.9972523, 0], [0.0052513, -0.0041175], [0.0074606, 0.0048125], [-0.0153783, -0.1968253], [0.0636871, -0.1408027], [0.3660976, -0.2937382]],
      gs48 = [[0.98879, 0], [0, 0], [-0.050909, 0], [0, 0], [0.075528, 0]],
      gs50 = [[0.9842990, 0], [0.0211642, 0.0037608], [-0.1036018, -0.0575102], [-0.0329095, -0.0320119], [0.0499471, 0.1223335], [0.0260460, 0.0899805], [0.0007388, -0.1435792], [0.0075848, -0.1334108], [-0.0216473, 0.0776645], [-0.0225161, 0.0853673]],
      miller$1 = [[0.9245, 0], [0, 0], [0.01943, 0]],
      lee = [[0.721316, 0], [0, 0], [-0.00881625, -0.00617325]];

  function modifiedStereographicAlaska() {
    return modifiedStereographic(alaska, [152, -64])
        .scale(1500)
        .center([-160.908, 62.4864])
        .clipAngle(25);
  }

  function modifiedStereographicGs48() {
    return modifiedStereographic(gs48, [95, -38])
        .scale(1000)
        .clipAngle(55)
        .center([-96.5563, 38.8675]);
  }

  function modifiedStereographicGs50() {
    return modifiedStereographic(gs50, [120, -45])
        .scale(359.513)
        .clipAngle(55)
        .center([-117.474, 53.0628]);
  }

  function modifiedStereographicMiller() {
    return modifiedStereographic(miller$1, [-20, -18])
        .scale(209.091)
        .center([20, 16.7214])
        .clipAngle(82);
  }

  function modifiedStereographicLee() {
    return modifiedStereographic(lee, [165, 10])
        .scale(250)
        .clipAngle(130)
        .center([-165, -10]);
  }

  function modifiedStereographic(coefficients, rotate) {
    var p = projection(modifiedStereographicRaw(coefficients)).rotate(rotate).clipAngle(90),
        r = rotation(rotate),
        center = p.center;

    delete p.rotate;

    p.center = function(_) {
      return arguments.length ? center(r(_)) : r.invert(center());
    };

    return p;
  }

  var sqrt6 = sqrt$1(6),
      sqrt7 = sqrt$1(7);

  function mtFlatPolarParabolicRaw(lambda, phi) {
    var theta = asin$1(7 * sin$1(phi) / (3 * sqrt6));
    return [
      sqrt6 * lambda * (2 * cos$1(2 * theta / 3) - 1) / sqrt7,
      9 * sin$1(theta / 3) / sqrt7
    ];
  }

  mtFlatPolarParabolicRaw.invert = function(x, y) {
    var theta = 3 * asin$1(y * sqrt7 / 9);
    return [
      x * sqrt7 / (sqrt6 * (2 * cos$1(2 * theta / 3) - 1)),
      asin$1(sin$1(theta) * 3 * sqrt6 / 7)
    ];
  };

  function mtFlatPolarParabolic() {
    return projection(mtFlatPolarParabolicRaw)
        .scale(164.859);
  }

  function mtFlatPolarQuarticRaw(lambda, phi) {
    var k = (1 + sqrt1_2) * sin$1(phi),
        theta = phi;
    for (var i = 0, delta; i < 25; i++) {
      theta -= delta = (sin$1(theta / 2) + sin$1(theta) - k) / (0.5 * cos$1(theta / 2) + cos$1(theta));
      if (abs$1(delta) < epsilon$1) break;
    }
    return [
      lambda * (1 + 2 * cos$1(theta) / cos$1(theta / 2)) / (3 * sqrt2),
      2 * sqrt$1(3) * sin$1(theta / 2) / sqrt$1(2 + sqrt2)
    ];
  }

  mtFlatPolarQuarticRaw.invert = function(x, y) {
    var sinTheta_2 = y * sqrt$1(2 + sqrt2) / (2 * sqrt$1(3)),
        theta = 2 * asin$1(sinTheta_2);
    return [
      3 * sqrt2 * x / (1 + 2 * cos$1(theta) / cos$1(theta / 2)),
      asin$1((sinTheta_2 + sin$1(theta)) / (1 + sqrt1_2))
    ];
  };

  function mtFlatPolarQuartic() {
    return projection(mtFlatPolarQuarticRaw)
        .scale(188.209);
  }

  function mtFlatPolarSinusoidalRaw(lambda, phi) {
    var A = sqrt$1(6 / (4 + pi$1)),
        k = (1 + pi$1 / 4) * sin$1(phi),
        theta = phi / 2;
    for (var i = 0, delta; i < 25; i++) {
      theta -= delta = (theta / 2 + sin$1(theta) - k) / (0.5 + cos$1(theta));
      if (abs$1(delta) < epsilon$1) break;
    }
    return [
      A * (0.5 + cos$1(theta)) * lambda / 1.5,
      A * theta
    ];
  }

  mtFlatPolarSinusoidalRaw.invert = function(x, y) {
    var A = sqrt$1(6 / (4 + pi$1)),
        theta = y / A;
    if (abs$1(abs$1(theta) - halfPi$1) < epsilon$1) theta = theta < 0 ? -halfPi$1 : halfPi$1;
    return [
      1.5 * x / (A * (0.5 + cos$1(theta))),
      asin$1((theta / 2 + sin$1(theta)) / (1 + pi$1 / 4))
    ];
  };

  function mtFlatPolarSinusoidal() {
    return projection(mtFlatPolarSinusoidalRaw)
        .scale(166.518);
  }

  function naturalEarth2Raw(lambda, phi) {
    var phi2 = phi * phi, phi4 = phi2 * phi2, phi6 = phi2 * phi4;
    return [
      lambda * (0.84719 - 0.13063 * phi2 + phi6 * phi6 * (-0.04515 + 0.05494 * phi2 - 0.02326 * phi4 + 0.00331 * phi6)),
      phi * (1.01183 + phi4 * phi4 * (-0.02625 + 0.01926 * phi2 - 0.00396 * phi4))
    ];
  }

  naturalEarth2Raw.invert = function(x, y) {
    var phi = y, i = 25, delta, phi2, phi4, phi6;
    do {
      phi2 = phi * phi; phi4 = phi2 * phi2;
      phi -= delta = ((phi * (1.01183 + phi4 * phi4 * (-0.02625 + 0.01926 * phi2 - 0.00396 * phi4))) - y) /
        (1.01183 + phi4 * phi4 * ((9 * -0.02625) + (11 * 0.01926) * phi2 + (13 * -0.00396) * phi4));
    } while (abs$1(delta) > epsilon2$1 && --i > 0);
    phi2 = phi * phi; phi4 = phi2 * phi2; phi6 = phi2 * phi4;
    return [
      x / (0.84719 - 0.13063 * phi2 + phi6 * phi6 * (-0.04515 + 0.05494 * phi2 - 0.02326 * phi4 + 0.00331 * phi6)),
      phi
    ];
  };

  function naturalEarth2() {
    return projection(naturalEarth2Raw)
        .scale(175.295);
  }

  function nellHammerRaw(lambda, phi) {
    return [
      lambda * (1 + cos$1(phi)) / 2,
      2 * (phi - tan$1(phi / 2))
    ];
  }

  nellHammerRaw.invert = function(x, y) {
    var p = y / 2;
    for (var i = 0, delta = Infinity; i < 10 && abs$1(delta) > epsilon$1; ++i) {
      var c = cos$1(y / 2);
      y -= delta = (y - tan$1(y / 2) - p) / (1 - 0.5 / (c * c));
    }
    return [
      2 * x / (1 + cos$1(y)),
      y
    ];
  };

  function nellHammer() {
    return projection(nellHammerRaw)
        .scale(152.63);
  }

  // Based on Java implementation by Bojan Savric.
  // https://github.com/OSUCartography/JMapProjLib/blob/master/src/com/jhlabs/map/proj/PattersonProjection.java

  var pattersonK1 = 1.0148,
      pattersonK2 = 0.23185,
      pattersonK3 = -0.14499,
      pattersonK4 = 0.02406,
      pattersonC1 = pattersonK1,
      pattersonC2 = 5 * pattersonK2,
      pattersonC3 = 7 * pattersonK3,
      pattersonC4 = 9 * pattersonK4,
      pattersonYmax = 1.790857183;

  function pattersonRaw(lambda, phi) {
    var phi2 = phi * phi;
    return [
      lambda,
      phi * (pattersonK1 + phi2 * phi2 * (pattersonK2 + phi2 * (pattersonK3 + pattersonK4 * phi2)))
    ];
  }

  pattersonRaw.invert = function(x, y) {
    if (y > pattersonYmax) y = pattersonYmax;
    else if (y < -pattersonYmax) y = -pattersonYmax;
    var yc = y, delta;

    do { // Newton-Raphson
      var y2 = yc * yc;
      yc -= delta = ((yc * (pattersonK1 + y2 * y2 * (pattersonK2 + y2 * (pattersonK3 + pattersonK4 * y2)))) - y) / (pattersonC1 + y2 * y2 * (pattersonC2 + y2 * (pattersonC3 + pattersonC4 * y2)));
    } while (abs$1(delta) > epsilon$1);

    return [x, yc];
  };

  function patterson() {
    return projection(pattersonRaw)
        .scale(139.319);
  }

  function polyconicRaw(lambda, phi) {
    if (abs$1(phi) < epsilon$1) return [lambda, 0];
    var tanPhi = tan$1(phi),
        k = lambda * sin$1(phi);
    return [
      sin$1(k) / tanPhi,
      phi + (1 - cos$1(k)) / tanPhi
    ];
  }

  polyconicRaw.invert = function(x, y) {
    if (abs$1(y) < epsilon$1) return [x, 0];
    var k = x * x + y * y,
        phi = y * 0.5,
        i = 10, delta;
    do {
      var tanPhi = tan$1(phi),
          secPhi = 1 / cos$1(phi),
          j = k - 2 * y * phi + phi * phi;
      phi -= delta = (tanPhi * j + 2 * (phi - y)) / (2 + j * secPhi * secPhi + 2 * (phi - y) * tanPhi);
    } while (abs$1(delta) > epsilon$1 && --i > 0);
    tanPhi = tan$1(phi);
    return [
      (abs$1(y) < abs$1(phi + 1 / tanPhi) ? asin$1(x * tanPhi) : sign$1(x) * (acos$1(abs$1(x * tanPhi)) + halfPi$1)) / sin$1(phi),
      phi
    ];
  };

  function polyconic() {
    return projection(polyconicRaw)
        .scale(103.74);
  }

  // Note: 6-element arrays are used to denote the 3x3 affine transform matrix:
  // [a, b, c,
  //  d, e, f,
  //  0, 0, 1] - this redundant row is left out.

  // Transform matrix for [a0, a1] -> [b0, b1].
  function matrix(a, b) {
    var u = subtract(a[1], a[0]),
        v = subtract(b[1], b[0]),
        phi = angle$2(u, v),
        s = length$2(u) / length$2(v);

    return multiply([
      1, 0, a[0][0],
      0, 1, a[0][1]
    ], multiply([
      s, 0, 0,
      0, s, 0
    ], multiply([
      cos$1(phi), sin$1(phi), 0,
      -sin$1(phi), cos$1(phi), 0
    ], [
      1, 0, -b[0][0],
      0, 1, -b[0][1]
    ])));
  }

  // Inverts a transform matrix.
  function inverse(m) {
    var k = 1 / (m[0] * m[4] - m[1] * m[3]);
    return [
      k * m[4], -k * m[1], k * (m[1] * m[5] - m[2] * m[4]),
      -k * m[3], k * m[0], k * (m[2] * m[3] - m[0] * m[5])
    ];
  }

  // Multiplies two 3x2 matrices.
  function multiply(a, b) {
    return [
      a[0] * b[0] + a[1] * b[3],
      a[0] * b[1] + a[1] * b[4],
      a[0] * b[2] + a[1] * b[5] + a[2],
      a[3] * b[0] + a[4] * b[3],
      a[3] * b[1] + a[4] * b[4],
      a[3] * b[2] + a[4] * b[5] + a[5]
    ];
  }

  // Subtracts 2D vectors.
  function subtract(a, b) {
    return [a[0] - b[0], a[1] - b[1]];
  }

  // Magnitude of a 2D vector.
  function length$2(v) {
    return sqrt$1(v[0] * v[0] + v[1] * v[1]);
  }

  // Angle between two 2D vectors.
  function angle$2(a, b) {
    return atan2$1(a[0] * b[1] - a[1] * b[0], a[0] * b[0] + a[1] * b[1]);
  }

  // Creates a polyhedral projection.
  //  * root: a spanning tree of polygon faces.  Nodes are automatically
  //    augmented with a transform matrix.
  //  * face: a function that returns the appropriate node for a given {lambda, phi}
  //    point (radians).
  //  * r: rotation angle for root face [deprecated by .angle()].
  function polyhedral(root, face, r) {

    recurse(root, {transform: null});

    function recurse(node, parent) {
      node.edges = faceEdges(node.face);
      // Find shared edge.
      if (parent.face) {
        var shared = node.shared = sharedEdge(node.face, parent.face),
            m = matrix(shared.map(parent.project), shared.map(node.project));
        node.transform = parent.transform ? multiply(parent.transform, m) : m;
        // Replace shared edge in parent edges array.
        var edges = parent.edges;
        for (var i = 0, n = edges.length; i < n; ++i) {
          if (pointEqual$2(shared[0], edges[i][1]) && pointEqual$2(shared[1], edges[i][0])) edges[i] = node;
          if (pointEqual$2(shared[0], edges[i][0]) && pointEqual$2(shared[1], edges[i][1])) edges[i] = node;
        }
        edges = node.edges;
        for (i = 0, n = edges.length; i < n; ++i) {
          if (pointEqual$2(shared[0], edges[i][0]) && pointEqual$2(shared[1], edges[i][1])) edges[i] = parent;
          if (pointEqual$2(shared[0], edges[i][1]) && pointEqual$2(shared[1], edges[i][0])) edges[i] = parent;
        }
      } else {
        node.transform = parent.transform;
      }
      if (node.children) {
        node.children.forEach(function(child) {
          recurse(child, node);
        });
      }
      return node;
    }

    function forward(lambda, phi) {
      var node = face(lambda, phi),
          point = node.project([lambda * degrees$1, phi * degrees$1]),
          t;
      if (t = node.transform) {
        return [
          t[0] * point[0] + t[1] * point[1] + t[2],
          -(t[3] * point[0] + t[4] * point[1] + t[5])
        ];
      }
      point[1] = -point[1];
      return point;
    }

    // Naive inverse!  A faster solution would use bounding boxes, or even a
    // polygonal quadtree.
    if (hasInverse(root)) forward.invert = function(x, y) {
      var coordinates = faceInvert(root, [x, -y]);
      return coordinates && (coordinates[0] *= radians$1, coordinates[1] *= radians$1, coordinates);
    };

    function faceInvert(node, coordinates) {
      var invert = node.project.invert,
          t = node.transform,
          point = coordinates;
      if (t) {
        t = inverse(t);
        point = [
          t[0] * point[0] + t[1] * point[1] + t[2],
          (t[3] * point[0] + t[4] * point[1] + t[5])
        ];
      }
      if (invert && node === faceDegrees(p = invert(point))) return p;
      var p,
          children = node.children;
      for (var i = 0, n = children && children.length; i < n; ++i) {
        if (p = faceInvert(children[i], coordinates)) return p;
      }
    }

    function faceDegrees(coordinates) {
      return face(coordinates[0] * radians$1, coordinates[1] * radians$1);
    }

    var proj = projection(forward),
        stream_ = proj.stream;

    proj.stream = function(stream) {
      var rotate = proj.rotate(),
          rotateStream = stream_(stream),
          sphereStream = (proj.rotate([0, 0]), stream_(stream));
      proj.rotate(rotate);
      rotateStream.sphere = function() {
        sphereStream.polygonStart();
        sphereStream.lineStart();
        outline(sphereStream, root);
        sphereStream.lineEnd();
        sphereStream.polygonEnd();
      };
      return rotateStream;
    };

    return proj.angle(r == null ? -30 : r * degrees$1);
  }

  function outline(stream, node, parent) {
    var point,
        edges = node.edges,
        n = edges.length,
        edge,
        multiPoint = {type: "MultiPoint", coordinates: node.face},
        notPoles = node.face.filter(function(d) { return abs$1(d[1]) !== 90; }),
        b = bounds({type: "MultiPoint", coordinates: notPoles}),
        inside = false,
        j = -1,
        dx = b[1][0] - b[0][0];
    // TODO
    var c = dx === 180 || dx === 360
        ? [(b[0][0] + b[1][0]) / 2, (b[0][1] + b[1][1]) / 2]
        : centroid(multiPoint);
    // First find the shared edge…
    if (parent) while (++j < n) {
      if (edges[j] === parent) break;
    }
    ++j;
    for (var i = 0; i < n; ++i) {
      edge = edges[(i + j) % n];
      if (Array.isArray(edge)) {
        if (!inside) {
          stream.point((point = interpolate(edge[0], c)(epsilon$1))[0], point[1]);
          inside = true;
        }
        stream.point((point = interpolate(edge[1], c)(epsilon$1))[0], point[1]);
      } else {
        inside = false;
        if (edge !== parent) outline(stream, edge, node);
      }
    }
  }

  // Tests equality of two spherical points.
  function pointEqual$2(a, b) {
    return a && b && a[0] === b[0] && a[1] === b[1];
  }

  // Finds a shared edge given two clockwise polygons.
  function sharedEdge(a, b) {
    var x, y, n = a.length, found = null;
    for (var i = 0; i < n; ++i) {
      x = a[i];
      for (var j = b.length; --j >= 0;) {
        y = b[j];
        if (x[0] === y[0] && x[1] === y[1]) {
          if (found) return [found, x];
          found = x;
        }
      }
    }
  }

  // Converts an array of n face vertices to an array of n + 1 edges.
  function faceEdges(face) {
    var n = face.length,
        edges = [];
    for (var a = face[n - 1], i = 0; i < n; ++i) edges.push([a, a = face[i]]);
    return edges;
  }

  function hasInverse(node) {
    return node.project.invert || node.children && node.children.some(hasInverse);
  }

  // TODO generate on-the-fly to avoid external modification.
  var octahedron = [
    [0, 90],
    [-90, 0], [0, 0], [90, 0], [180, 0],
    [0, -90]
  ];

  var octahedron$1 = [
    [0, 2, 1],
    [0, 3, 2],
    [5, 1, 2],
    [5, 2, 3],
    [0, 1, 4],
    [0, 4, 3],
    [5, 4, 1],
    [5, 3, 4]
  ].map(function(face) {
    return face.map(function(i) {
      return octahedron[i];
    });
  });

  function butterfly(faceProjection) {

    faceProjection = faceProjection || function(face) {
      var c = centroid({type: "MultiPoint", coordinates: face});
      return gnomonic().scale(1).translate([0, 0]).rotate([-c[0], -c[1]]);
    };

    var faces = octahedron$1.map(function(face) {
      return {face: face, project: faceProjection(face)};
    });

    [-1, 0, 0, 1, 0, 1, 4, 5].forEach(function(d, i) {
      var node = faces[d];
      node && (node.children || (node.children = [])).push(faces[i]);
    });

    return polyhedral(faces[0], function(lambda, phi) {
          return faces[lambda < -pi$1 / 2 ? phi < 0 ? 6 : 4
              : lambda < 0 ? phi < 0 ? 2 : 0
              : lambda < pi$1 / 2 ? phi < 0 ? 3 : 1
              : phi < 0 ? 7 : 5];
        })
        .angle(-30)
        .scale(101.858)
        .center([0, 45]);
  }

  var kx = 2 / sqrt$1(3);

  function collignonK(a, b) {
    var p = collignonRaw(a, b);
    return [p[0] * kx, p[1]];
  }

  collignonK.invert = function(x,y) {
    return collignonRaw.invert(x / kx, y);
  };

  function collignon$1(faceProjection) {

    faceProjection = faceProjection || function(face) {
      var c = centroid({type: "MultiPoint", coordinates: face});
      return projection(collignonK).translate([0, 0]).scale(1).rotate(c[1] > 0 ? [-c[0], 0] : [180 - c[0], 180]);
    };

    var faces = octahedron$1.map(function(face) {
      return {face: face, project: faceProjection(face)};
    });

    [-1, 0, 0, 1, 0, 1, 4, 5].forEach(function(d, i) {
      var node = faces[d];
      node && (node.children || (node.children = [])).push(faces[i]);
    });

    return polyhedral(faces[0], function(lambda, phi) {
          return faces[lambda < -pi$1 / 2 ? phi < 0 ? 6 : 4
              : lambda < 0 ? phi < 0 ? 2 : 0
              : lambda < pi$1 / 2 ? phi < 0 ? 3 : 1
              : phi < 0 ? 7 : 5];
        })
        .angle(-30)
        .scale(121.906)
        .center([0, 48.5904]);
  }

  function waterman(faceProjection) {

    faceProjection = faceProjection || function(face) {
      var c = face.length === 6 ? centroid({type: "MultiPoint", coordinates: face}) : face[0];
      return gnomonic().scale(1).translate([0, 0]).rotate([-c[0], -c[1]]);
    };

    var w5 = octahedron$1.map(function(face) {
      var xyz = face.map(cartesian$1),
          n = xyz.length,
          a = xyz[n - 1],
          b,
          hexagon = [];
      for (var i = 0; i < n; ++i) {
        b = xyz[i];
        hexagon.push(spherical$1([
          a[0] * 0.9486832980505138 + b[0] * 0.31622776601683794,
          a[1] * 0.9486832980505138 + b[1] * 0.31622776601683794,
          a[2] * 0.9486832980505138 + b[2] * 0.31622776601683794
        ]), spherical$1([
          b[0] * 0.9486832980505138 + a[0] * 0.31622776601683794,
          b[1] * 0.9486832980505138 + a[1] * 0.31622776601683794,
          b[2] * 0.9486832980505138 + a[2] * 0.31622776601683794
        ]));
        a = b;
      }
      return hexagon;
    });

    var cornerNormals = [];

    var parents = [-1, 0, 0, 1, 0, 1, 4, 5];

    w5.forEach(function(hexagon, j) {
      var face = octahedron$1[j],
          n = face.length,
          normals = cornerNormals[j] = [];
      for (var i = 0; i < n; ++i) {
        w5.push([
          face[i],
          hexagon[(i * 2 + 2) % (2 * n)],
          hexagon[(i * 2 + 1) % (2 * n)]
        ]);
        parents.push(j);
        normals.push(cross$1(
          cartesian$1(hexagon[(i * 2 + 2) % (2 * n)]),
          cartesian$1(hexagon[(i * 2 + 1) % (2 * n)])
        ));
      }
    });

    var faces = w5.map(function(face) {
      return {
        project: faceProjection(face),
        face: face
      };
    });

    parents.forEach(function(d, i) {
      var parent = faces[d];
      parent && (parent.children || (parent.children = [])).push(faces[i]);
    });

    function face(lambda, phi) {
      var cosphi = cos$1(phi),
          p = [cosphi * cos$1(lambda), cosphi * sin$1(lambda), sin$1(phi)];

      var hexagon = lambda < -pi$1 / 2 ? phi < 0 ? 6 : 4
          : lambda < 0 ? phi < 0 ? 2 : 0
          : lambda < pi$1 / 2 ? phi < 0 ? 3 : 1
          : phi < 0 ? 7 : 5;

      var n = cornerNormals[hexagon];

      return faces[dot(n[0], p) < 0 ? 8 + 3 * hexagon
          : dot(n[1], p) < 0 ? 8 + 3 * hexagon + 1
          : dot(n[2], p) < 0 ? 8 + 3 * hexagon + 2
          : hexagon];
    }

    return polyhedral(faces[0], face)
        .angle(-30)
        .scale(110.625)
        .center([0,45]);
  }

  function dot(a, b) {
    for (var i = 0, n = a.length, s = 0; i < n; ++i) s += a[i] * b[i];
    return s;
  }

  function cross$1(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]
    ];
  }

  // Converts 3D Cartesian to spherical coordinates (degrees).
  function spherical$1(cartesian) {
    return [
      atan2$1(cartesian[1], cartesian[0]) * degrees$1,
      asin$1(max$1(-1, min$1(1, cartesian[2]))) * degrees$1
    ];
  }

  // Converts spherical coordinates (degrees) to 3D Cartesian.
  function cartesian$1(coordinates) {
    var lambda = coordinates[0] * radians$1,
        phi = coordinates[1] * radians$1,
        cosphi = cos$1(phi);
    return [
      cosphi * cos$1(lambda),
      cosphi * sin$1(lambda),
      sin$1(phi)
    ];
  }

  function noop$1() {}

  function clockwise(ring) {
    if ((n = ring.length) < 4) return false;
    var i = 0,
        n,
        area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
    while (++i < n) area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
    return area <= 0;
  }

  function contains$1(ring, point) {
    var x = point[0],
        y = point[1],
        contains = false;
    for (var i = 0, n = ring.length, j = n - 1; i < n; j = i++) {
      var pi = ring[i], xi = pi[0], yi = pi[1],
          pj = ring[j], xj = pj[0], yj = pj[1];
      if (((yi > y) ^ (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) contains = !contains;
    }
    return contains;
  }

  function index$1(object, projection$$1) {
    var stream = projection$$1.stream, project;
    if (!stream) throw new Error("invalid projection");
    switch (object && object.type) {
      case "Feature": project = projectFeature; break;
      case "FeatureCollection": project = projectFeatureCollection; break;
      default: project = projectGeometry; break;
    }
    return project(object, stream);
  }

  function projectFeatureCollection(o, stream) {
    return {
      type: "FeatureCollection",
      features: o.features.map(function(f) {
        return projectFeature(f, stream);
      })
    };
  }

  function projectFeature(o, stream) {
    return {
      type: "Feature",
      id: o.id,
      properties: o.properties,
      geometry: projectGeometry(o.geometry, stream)
    };
  }

  function projectGeometryCollection(o, stream) {
    return {
      type: "GeometryCollection",
      geometries: o.geometries.map(function(o) {
        return projectGeometry(o, stream);
      })
    };
  }

  function projectGeometry(o, stream) {
    if (!o) return null;
    if (o.type === "GeometryCollection") return projectGeometryCollection(o, stream);
    var sink;
    switch (o.type) {
      case "Point": sink = sinkPoint; break;
      case "MultiPoint": sink = sinkPoint; break;
      case "LineString": sink = sinkLine; break;
      case "MultiLineString": sink = sinkLine; break;
      case "Polygon": sink = sinkPolygon; break;
      case "MultiPolygon": sink = sinkPolygon; break;
      case "Sphere": sink = sinkPolygon; break;
      default: return null;
    }
    geoStream(o, stream(sink));
    return sink.result();
  }

  var points = [],
      lines = [];

  var sinkPoint = {
    point: function(x, y) {
      points.push([x, y]);
    },
    result: function() {
      var result = !points.length ? null
          : points.length < 2 ? {type: "Point", coordinates: points[0]}
          : {type: "MultiPoint", coordinates: points};
      points = [];
      return result;
    }
  };

  var sinkLine = {
    lineStart: noop$1,
    point: function(x, y) {
      points.push([x, y]);
    },
    lineEnd: function() {
      if (points.length) lines.push(points), points = [];
    },
    result: function() {
      var result = !lines.length ? null
          : lines.length < 2 ? {type: "LineString", coordinates: lines[0]}
          : {type: "MultiLineString", coordinates: lines};
      lines = [];
      return result;
    }
  };

  var sinkPolygon = {
    polygonStart: noop$1,
    lineStart: noop$1,
    point: function(x, y) {
      points.push([x, y]);
    },
    lineEnd: function() {
      var n = points.length;
      if (n) {
        do points.push(points[0].slice()); while (++n < 4);
        lines.push(points), points = [];
      }
    },
    polygonEnd: noop$1,
    result: function() {
      if (!lines.length) return null;
      var polygons = [],
          holes = [];

      // https://github.com/d3/d3/issues/1558
      lines.forEach(function(ring) {
        if (clockwise(ring)) polygons.push([ring]);
        else holes.push(ring);
      });

      holes.forEach(function(hole) {
        var point = hole[0];
        polygons.some(function(polygon) {
          if (contains$1(polygon[0], point)) {
            polygon.push(hole);
            return true;
          }
        }) || polygons.push([hole]);
      });

      lines = [];

      return !polygons.length ? null
          : polygons.length > 1 ? {type: "MultiPolygon", coordinates: polygons}
          : {type: "Polygon", coordinates: polygons[0]};
    }
  };

  function quincuncial(project) {
    var dx = project(halfPi$1, 0)[0] - project(-halfPi$1, 0)[0];

    function projectQuincuncial(lambda, phi) {
      var t = abs$1(lambda) < halfPi$1,
          p = project(t ? lambda : lambda > 0 ? lambda - pi$1 : lambda + pi$1, phi),
          x = (p[0] - p[1]) * sqrt1_2,
          y = (p[0] + p[1]) * sqrt1_2;
      if (t) return [x, y];
      var d = dx * sqrt1_2,
          s = x > 0 ^ y > 0 ? -1 : 1;
      return [s * x - sign$1(y) * d, s * y - sign$1(x) * d];
    }

    if (project.invert) projectQuincuncial.invert = function(x0, y0) {
      var x = (x0 + y0) * sqrt1_2,
          y = (y0 - x0) * sqrt1_2,
          t = abs$1(x) < 0.5 * dx && abs$1(y) < 0.5 * dx;

      if (!t) {
        var d = dx * sqrt1_2,
            s = x > 0 ^ y > 0 ? -1 : 1,
            x1 = -s * x0 + (y > 0 ? 1 : -1) * d,
            y1 = -s * y0 + (x > 0 ? 1 : -1) * d;
        x = (-x1 - y1) * sqrt1_2;
        y = (x1 - y1) * sqrt1_2;
      }

      var p = project.invert(x, y);
      if (!t) p[0] += x > 0 ? pi$1 : -pi$1;
      return p;
    };

    return projection(projectQuincuncial)
        .rotate([-90, -90, 45])
        .clipAngle(180 - 1e-3);
  }

  function gringorten$1() {
    return quincuncial(gringortenRaw)
        .scale(176.423);
  }

  function peirce() {
    return quincuncial(guyouRaw)
        .scale(111.48);
  }

  function quantize(input, digits) {
    if (!(0 <= (digits = +digits) && digits <= 20)) throw new Error("invalid digits");

    function quantizePoint(input) {
      var n = input.length, i = 2, output = new Array(n);
      output[0] = +input[0].toFixed(digits);
      output[1] = +input[1].toFixed(digits);
      while (i < n) output[i] = input[i], ++i;
      return output;
    }

    function quantizePoints(input) {
      return input.map(quantizePoint);
    }

    function quantizePolygon(input) {
      return input.map(quantizePoints);
    }

    function quantizeGeometry(input) {
      if (input == null) return input;
      var output;
      switch (input.type) {
        case "GeometryCollection": output = {type: "GeometryCollection", geometries: input.geometries.map(quantizeGeometry)}; break;
        case "Point": output = {type: "Point", coordinates: quantizePoint(input.coordinates)}; break;
        case "MultiPoint": case "LineString": output = {type: input.type, coordinates: quantizePoints(input.coordinates)}; break;
        case "MultiLineString": case "Polygon": output = {type: input.type, coordinates: quantizePolygon(input.coordinates)}; break;
        case "MultiPolygon": output = {type: "MultiPolygon", coordinates: input.coordinates.map(quantizePolygon)}; break;
        default: return input;
      }
      if (input.bbox != null) output.bbox = input.bbox;
      return output;
    }

    function quantizeFeature(input) {
      var output = {type: "Feature", properties: input.properties, geometry: quantizeGeometry(input.geometry)};
      if (input.id != null) output.id = input.id;
      if (input.bbox != null) output.bbox = input.bbox;
      return output;
    }

    if (input != null) switch (input.type) {
      case "Feature": return quantizeFeature(input);
      case "FeatureCollection": {
        var output = {type: "FeatureCollection", features: input.features.map(quantizeFeature)};
        if (input.bbox != null) output.bbox = input.bbox;
        return output;
      }
      default: return quantizeGeometry(input);
    }

    return input;
  }

  function rectangularPolyconicRaw(phi0) {
    var sinPhi0 = sin$1(phi0);

    function forward(lambda, phi) {
      var A = sinPhi0 ? tan$1(lambda * sinPhi0 / 2) / sinPhi0 : lambda / 2;
      if (!phi) return [2 * A, -phi0];
      var E = 2 * atan$1(A * sin$1(phi)),
          cotPhi = 1 / tan$1(phi);
      return [
        sin$1(E) * cotPhi,
        phi + (1 - cos$1(E)) * cotPhi - phi0
      ];
    }

    // TODO return null for points outside outline.
    forward.invert = function(x, y) {
      if (abs$1(y += phi0) < epsilon$1) return [sinPhi0 ? 2 * atan$1(sinPhi0 * x / 2) / sinPhi0 : x, 0];
      var k = x * x + y * y,
          phi = 0,
          i = 10, delta;
      do {
        var tanPhi = tan$1(phi),
            secPhi = 1 / cos$1(phi),
            j = k - 2 * y * phi + phi * phi;
        phi -= delta = (tanPhi * j + 2 * (phi - y)) / (2 + j * secPhi * secPhi + 2 * (phi - y) * tanPhi);
      } while (abs$1(delta) > epsilon$1 && --i > 0);
      var E = x * (tanPhi = tan$1(phi)),
          A = tan$1(abs$1(y) < abs$1(phi + 1 / tanPhi) ? asin$1(E) * 0.5 : acos$1(E) * 0.5 + pi$1 / 4) / sin$1(phi);
      return [
        sinPhi0 ? 2 * atan$1(sinPhi0 * A) / sinPhi0 : 2 * A,
        phi
      ];
    };

    return forward;
  }

  function rectangularPolyconic() {
    return parallel1(rectangularPolyconicRaw)
        .scale(131.215);
  }

  var K = [
    [0.9986, -0.062],
    [1.0000, 0.0000],
    [0.9986, 0.0620],
    [0.9954, 0.1240],
    [0.9900, 0.1860],
    [0.9822, 0.2480],
    [0.9730, 0.3100],
    [0.9600, 0.3720],
    [0.9427, 0.4340],
    [0.9216, 0.4958],
    [0.8962, 0.5571],
    [0.8679, 0.6176],
    [0.8350, 0.6769],
    [0.7986, 0.7346],
    [0.7597, 0.7903],
    [0.7186, 0.8435],
    [0.6732, 0.8936],
    [0.6213, 0.9394],
    [0.5722, 0.9761],
    [0.5322, 1.0000]
  ];

  K.forEach(function(d) {
    d[1] *= 1.0144;
  });

  function robinsonRaw(lambda, phi) {
    var i = min$1(18, abs$1(phi) * 36 / pi$1),
        i0 = floor$1(i),
        di = i - i0,
        ax = (k = K[i0])[0],
        ay = k[1],
        bx = (k = K[++i0])[0],
        by = k[1],
        cx = (k = K[min$1(19, ++i0)])[0],
        cy = k[1],
        k;
    return [
      lambda * (bx + di * (cx - ax) / 2 + di * di * (cx - 2 * bx + ax) / 2),
      (phi > 0 ? halfPi$1 : -halfPi$1) * (by + di * (cy - ay) / 2 + di * di * (cy - 2 * by + ay) / 2)
    ];
  }

  robinsonRaw.invert = function(x, y) {
    var yy = y / halfPi$1,
        phi = yy * 90,
        i = min$1(18, abs$1(phi / 5)),
        i0 = max$1(0, floor$1(i));
    do {
      var ay = K[i0][1],
          by = K[i0 + 1][1],
          cy = K[min$1(19, i0 + 2)][1],
          u = cy - ay,
          v = cy - 2 * by + ay,
          t = 2 * (abs$1(yy) - by) / u,
          c = v / u,
          di = t * (1 - c * t * (1 - 2 * c * t));
      if (di >= 0 || i0 === 1) {
        phi = (y >= 0 ? 5 : -5) * (di + i);
        var j = 50, delta;
        do {
          i = min$1(18, abs$1(phi) / 5);
          i0 = floor$1(i);
          di = i - i0;
          ay = K[i0][1];
          by = K[i0 + 1][1];
          cy = K[min$1(19, i0 + 2)][1];
          phi -= (delta = (y >= 0 ? halfPi$1 : -halfPi$1) * (by + di * (cy - ay) / 2 + di * di * (cy - 2 * by + ay) / 2) - y) * degrees$1;
        } while (abs$1(delta) > epsilon2$1 && --j > 0);
        break;
      }
    } while (--i0 >= 0);
    var ax = K[i0][0],
        bx = K[i0 + 1][0],
        cx = K[min$1(19, i0 + 2)][0];
    return [
      x / (bx + di * (cx - ax) / 2 + di * di * (cx - 2 * bx + ax) / 2),
      phi * radians$1
    ];
  };

  function robinson() {
    return projection(robinsonRaw)
        .scale(152.63);
  }

  function satelliteVerticalRaw(P) {
    function forward(lambda, phi) {
      var cosPhi = cos$1(phi),
          k = (P - 1) / (P - cosPhi * cos$1(lambda));
      return [
        k * cosPhi * sin$1(lambda),
        k * sin$1(phi)
      ];
    }

    forward.invert = function(x, y) {
      var rho2 = x * x + y * y,
          rho = sqrt$1(rho2),
          sinc = (P - sqrt$1(1 - rho2 * (P + 1) / (P - 1))) / ((P - 1) / rho + rho / (P - 1));
      return [
        atan2$1(x * sinc, rho * sqrt$1(1 - sinc * sinc)),
        rho ? asin$1(y * sinc / rho) : 0
      ];
    };

    return forward;
  }

  function satelliteRaw(P, omega) {
    var vertical = satelliteVerticalRaw(P);
    if (!omega) return vertical;
    var cosOmega = cos$1(omega),
        sinOmega = sin$1(omega);

    function forward(lambda, phi) {
      var coordinates = vertical(lambda, phi),
          y = coordinates[1],
          A = y * sinOmega / (P - 1) + cosOmega;
      return [
        coordinates[0] * cosOmega / A,
        y / A
      ];
    }

    forward.invert = function(x, y) {
      var k = (P - 1) / (P - 1 - y * sinOmega);
      return vertical.invert(k * x, k * y * cosOmega);
    };

    return forward;
  }

  function satellite() {
    var distance$$1 = 2,
        omega = 0,
        m = projectionMutator(satelliteRaw),
        p = m(distance$$1, omega);

    // As a multiple of radius.
    p.distance = function(_) {
      if (!arguments.length) return distance$$1;
      return m(distance$$1 = +_, omega);
    };

    p.tilt = function(_) {
      if (!arguments.length) return omega * degrees$1;
      return m(distance$$1, omega = _ * radians$1);
    };

    return p
        .scale(432.147)
        .clipAngle(acos$1(1 / distance$$1) * degrees$1 - 1e-6);
  }

  var epsilon$2 = 1e-4,
      epsilonInverse = 1e4,
      x0$5 = -180, x0e = x0$5 + epsilon$2,
      x1$1 = 180, x1e = x1$1 - epsilon$2,
      y0$5 = -90, y0e = y0$5 + epsilon$2,
      y1$1 = 90, y1e = y1$1 - epsilon$2;

  function nonempty(coordinates) {
    return coordinates.length > 0;
  }

  function quantize$1(x) {
    return Math.floor(x * epsilonInverse) / epsilonInverse;
  }

  function normalizePoint(y) {
    return y === y0$5 || y === y1$1 ? [0, y] : [x0$5, quantize$1(y)]; // pole or antimeridian?
  }

  function clampPoint(p) {
    var x = p[0], y = p[1], clamped = false;
    if (x <= x0e) x = x0$5, clamped = true;
    else if (x >= x1e) x = x1$1, clamped = true;
    if (y <= y0e) y = y0$5, clamped = true;
    else if (y >= y1e) y = y1$1, clamped = true;
    return clamped ? [x, y] : p;
  }

  function clampPoints(points) {
    return points.map(clampPoint);
  }

  // For each ring, detect where it crosses the antimeridian or pole.
  function extractFragments(rings, polygon, fragments) {
    for (var j = 0, m = rings.length; j < m; ++j) {
      var ring = rings[j].slice();

      // By default, assume that this ring doesn’t need any stitching.
      fragments.push({index: -1, polygon: polygon, ring: ring});

      for (var i = 0, n = ring.length; i < n; ++i) {
        var point = ring[i],
            x = point[0],
            y = point[1];

        // If this is an antimeridian or polar point…
        if (x <= x0e || x >= x1e || y <= y0e || y >= y1e) {
          ring[i] = clampPoint(point);

          // Advance through any antimeridian or polar points…
          for (var k = i + 1; k < n; ++k) {
            var pointk = ring[k],
                xk = pointk[0],
                yk = pointk[1];
            if (xk > x0e && xk < x1e && yk > y0e && yk < y1e) break;
          }

          // If this was just a single antimeridian or polar point,
          // we don’t need to cut this ring into a fragment;
          // we can just leave it as-is.
          if (k === i + 1) continue;

          // Otherwise, if this is not the first point in the ring,
          // cut the current fragment so that it ends at the current point.
          // The current point is also normalized for later joining.
          if (i) {
            var fragmentBefore = {index: -1, polygon: polygon, ring: ring.slice(0, i + 1)};
            fragmentBefore.ring[fragmentBefore.ring.length - 1] = normalizePoint(y);
            fragments[fragments.length - 1] = fragmentBefore;
          }

          // If the ring started with an antimeridian fragment,
          // we can ignore that fragment entirely.
          else fragments.pop();

          // If the remainder of the ring is an antimeridian fragment,
          // move on to the next ring.
          if (k >= n) break;

          // Otherwise, add the remaining ring fragment and continue.
          fragments.push({index: -1, polygon: polygon, ring: ring = ring.slice(k - 1)});
          ring[0] = normalizePoint(ring[0][1]);
          i = -1;
          n = ring.length;
        }
      }
    }
  }

  // Now stitch the fragments back together into rings.
  function stitchFragments(fragments) {
    var i, n = fragments.length;

    // To connect the fragments start-to-end, create a simple index by end.
    var fragmentByStart = {},
        fragmentByEnd = {},
        fragment,
        start,
        startFragment,
        end,
        endFragment;

    // For each fragment…
    for (i = 0; i < n; ++i) {
      fragment = fragments[i];
      start = fragment.ring[0];
      end = fragment.ring[fragment.ring.length - 1];

      // If this fragment is closed, add it as a standalone ring.
      if (start[0] === end[0] && start[1] === end[1]) {
        fragment.polygon.push(fragment.ring);
        fragments[i] = null;
        continue;
      }

      fragment.index = i;
      fragmentByStart[start] = fragmentByEnd[end] = fragment;
    }

    // For each open fragment…
    for (i = 0; i < n; ++i) {
      fragment = fragments[i];
      if (fragment) {
        start = fragment.ring[0];
        end = fragment.ring[fragment.ring.length - 1];
        startFragment = fragmentByEnd[start];
        endFragment = fragmentByStart[end];

        delete fragmentByStart[start];
        delete fragmentByEnd[end];

        // If this fragment is closed, add it as a standalone ring.
        if (start[0] === end[0] && start[1] === end[1]) {
          fragment.polygon.push(fragment.ring);
          continue;
        }

        if (startFragment) {
          delete fragmentByEnd[start];
          delete fragmentByStart[startFragment.ring[0]];
          startFragment.ring.pop(); // drop the shared coordinate
          fragments[startFragment.index] = null;
          fragment = {index: -1, polygon: startFragment.polygon, ring: startFragment.ring.concat(fragment.ring)};

          if (startFragment === endFragment) {
            // Connect both ends to this single fragment to create a ring.
            fragment.polygon.push(fragment.ring);
          } else {
            fragment.index = n++;
            fragments.push(fragmentByStart[fragment.ring[0]] = fragmentByEnd[fragment.ring[fragment.ring.length - 1]] = fragment);
          }
        } else if (endFragment) {
          delete fragmentByStart[end];
          delete fragmentByEnd[endFragment.ring[endFragment.ring.length - 1]];
          fragment.ring.pop(); // drop the shared coordinate
          fragment = {index: n++, polygon: endFragment.polygon, ring: fragment.ring.concat(endFragment.ring)};
          fragments[endFragment.index] = null;
          fragments.push(fragmentByStart[fragment.ring[0]] = fragmentByEnd[fragment.ring[fragment.ring.length - 1]] = fragment);
        } else {
          fragment.ring.push(fragment.ring[0]); // close ring
          fragment.polygon.push(fragment.ring);
        }
      }
    }
  }

  function stitchFeature(input) {
    var output = {type: "Feature", geometry: stitchGeometry(input.geometry)};
    if (input.id != null) output.id = input.id;
    if (input.bbox != null) output.bbox = input.bbox;
    if (input.properties != null) output.properties = input.properties;
    return output;
  }

  function stitchGeometry(input) {
    if (input == null) return input;
    var output, fragments, i, n;
    switch (input.type) {
      case "GeometryCollection": output = {type: "GeometryCollection", geometries: input.geometries.map(stitchGeometry)}; break;
      case "Point": output = {type: "Point", coordinates: clampPoint(input.coordinates)}; break;
      case "MultiPoint": case "LineString": output = {type: input.type, coordinates: clampPoints(input.coordinates)}; break;
      case "MultiLineString": output = {type: "MultiLineString", coordinates: input.coordinates.map(clampPoints)}; break;
      case "Polygon": {
        var polygon = [];
        extractFragments(input.coordinates, polygon, fragments = []);
        stitchFragments(fragments);
        output = {type: "Polygon", coordinates: polygon};
        break;
      }
      case "MultiPolygon": {
        fragments = [], i = -1, n = input.coordinates.length;
        var polygons = new Array(n);
        while (++i < n) extractFragments(input.coordinates[i], polygons[i] = [], fragments);
        stitchFragments(fragments);
        output = {type: "MultiPolygon", coordinates: polygons.filter(nonempty)};
        break;
      }
      default: return input;
    }
    if (input.bbox != null) output.bbox = input.bbox;
    return output;
  }

  function stitch(input) {
    if (input == null) return input;
    switch (input.type) {
      case "Feature": return stitchFeature(input);
      case "FeatureCollection": {
        var output = {type: "FeatureCollection", features: input.features.map(stitchFeature)};
        if (input.bbox != null) output.bbox = input.bbox;
        return output;
      }
      default: return stitchGeometry(input);
    }
  }

  function timesRaw(lambda, phi) {
    var t = tan$1(phi / 2),
        s = sin$1(quarterPi$1 * t);
    return [
      lambda * (0.74482 - 0.34588 * s * s),
      1.70711 * t
    ];
  }

  timesRaw.invert = function(x, y) {
    var t = y / 1.70711,
        s = sin$1(quarterPi$1 * t);
    return [
      x / (0.74482 - 0.34588 * s * s),
      2 * atan$1(t)
    ];
  };

  function times() {
    return projection(timesRaw)
        .scale(146.153);
  }

  // Compute the origin as the midpoint of the two reference points.
  // Rotate one of the reference points by the origin.
  // Apply the spherical law of sines to compute gamma rotation.
  function twoPoint(raw, p0, p1) {
    var i = interpolate(p0, p1),
        o = i(0.5),
        a = rotation([-o[0], -o[1]])(p0),
        b = i.distance / 2,
        y = -asin$1(sin$1(a[1] * radians$1) / sin$1(b)),
        R = [-o[0], -o[1], -(a[0] > 0 ? pi$1 - y : y) * degrees$1],
        p = projection(raw(b)).rotate(R),
        r = rotation(R),
        center = p.center;

    delete p.rotate;

    p.center = function(_) {
      return arguments.length ? center(r(_)) : r.invert(center());
    };

    return p
        .clipAngle(90);
  }

  function twoPointAzimuthalRaw(d) {
    var cosd = cos$1(d);

    function forward(lambda, phi) {
      var coordinates = gnomonicRaw(lambda, phi);
      coordinates[0] *= cosd;
      return coordinates;
    }

    forward.invert = function(x, y) {
      return gnomonicRaw.invert(x / cosd, y);
    };

    return forward;
  }

  function twoPointAzimuthalUsa() {
    return twoPointAzimuthal([-158, 21.5], [-77, 39])
        .clipAngle(60)
        .scale(400);
  }

  function twoPointAzimuthal(p0, p1) {
    return twoPoint(twoPointAzimuthalRaw, p0, p1);
  }

  // TODO clip to ellipse
  function twoPointEquidistantRaw(z0) {
    if (!(z0 *= 2)) return azimuthalEquidistantRaw;
    var lambdaa = -z0 / 2,
        lambdab = -lambdaa,
        z02 = z0 * z0,
        tanLambda0 = tan$1(lambdab),
        S = 0.5 / sin$1(lambdab);

    function forward(lambda, phi) {
      var za = acos$1(cos$1(phi) * cos$1(lambda - lambdaa)),
          zb = acos$1(cos$1(phi) * cos$1(lambda - lambdab)),
          ys = phi < 0 ? -1 : 1;
      za *= za, zb *= zb;
      return [
        (za - zb) / (2 * z0),
        ys * sqrt$1(4 * z02 * zb - (z02 - za + zb) * (z02 - za + zb)) / (2 * z0)
      ];
    }

    forward.invert = function(x, y) {
      var y2 = y * y,
          cosza = cos$1(sqrt$1(y2 + (t = x + lambdaa) * t)),
          coszb = cos$1(sqrt$1(y2 + (t = x + lambdab) * t)),
          t,
          d;
      return [
        atan2$1(d = cosza - coszb, t = (cosza + coszb) * tanLambda0),
        (y < 0 ? -1 : 1) * acos$1(sqrt$1(t * t + d * d) * S)
      ];
    };

    return forward;
  }

  function twoPointEquidistantUsa() {
    return twoPointEquidistant([-158, 21.5], [-77, 39])
        .clipAngle(130)
        .scale(122.571);
  }

  function twoPointEquidistant(p0, p1) {
    return twoPoint(twoPointEquidistantRaw, p0, p1);
  }

  function vanDerGrintenRaw(lambda, phi) {
    if (abs$1(phi) < epsilon$1) return [lambda, 0];
    var sinTheta = abs$1(phi / halfPi$1),
        theta = asin$1(sinTheta);
    if (abs$1(lambda) < epsilon$1 || abs$1(abs$1(phi) - halfPi$1) < epsilon$1) return [0, sign$1(phi) * pi$1 * tan$1(theta / 2)];
    var cosTheta = cos$1(theta),
        A = abs$1(pi$1 / lambda - lambda / pi$1) / 2,
        A2 = A * A,
        G = cosTheta / (sinTheta + cosTheta - 1),
        P = G * (2 / sinTheta - 1),
        P2 = P * P,
        P2_A2 = P2 + A2,
        G_P2 = G - P2,
        Q = A2 + G;
    return [
      sign$1(lambda) * pi$1 * (A * G_P2 + sqrt$1(A2 * G_P2 * G_P2 - P2_A2 * (G * G - P2))) / P2_A2,
      sign$1(phi) * pi$1 * (P * Q - A * sqrt$1((A2 + 1) * P2_A2 - Q * Q)) / P2_A2
    ];
  }

  vanDerGrintenRaw.invert = function(x, y) {
    if (abs$1(y) < epsilon$1) return [x, 0];
    if (abs$1(x) < epsilon$1) return [0, halfPi$1 * sin$1(2 * atan$1(y / pi$1))];
    var x2 = (x /= pi$1) * x,
        y2 = (y /= pi$1) * y,
        x2_y2 = x2 + y2,
        z = x2_y2 * x2_y2,
        c1 = -abs$1(y) * (1 + x2_y2),
        c2 = c1 - 2 * y2 + x2,
        c3 = -2 * c1 + 1 + 2 * y2 + z,
        d = y2 / c3 + (2 * c2 * c2 * c2 / (c3 * c3 * c3) - 9 * c1 * c2 / (c3 * c3)) / 27,
        a1 = (c1 - c2 * c2 / (3 * c3)) / c3,
        m1 = 2 * sqrt$1(-a1 / 3),
        theta1 = acos$1(3 * d / (a1 * m1)) / 3;
    return [
      pi$1 * (x2_y2 - 1 + sqrt$1(1 + 2 * (x2 - y2) + z)) / (2 * x),
      sign$1(y) * pi$1 * (-m1 * cos$1(theta1 + pi$1 / 3) - c2 / (3 * c3))
    ];
  };

  function vanDerGrinten() {
    return projection(vanDerGrintenRaw)
        .scale(79.4183);
  }

  function vanDerGrinten2Raw(lambda, phi) {
    if (abs$1(phi) < epsilon$1) return [lambda, 0];
    var sinTheta = abs$1(phi / halfPi$1),
        theta = asin$1(sinTheta);
    if (abs$1(lambda) < epsilon$1 || abs$1(abs$1(phi) - halfPi$1) < epsilon$1) return [0, sign$1(phi) * pi$1 * tan$1(theta / 2)];
    var cosTheta = cos$1(theta),
        A = abs$1(pi$1 / lambda - lambda / pi$1) / 2,
        A2 = A * A,
        x1 = cosTheta * (sqrt$1(1 + A2) - A * cosTheta) / (1 + A2 * sinTheta * sinTheta);
    return [
      sign$1(lambda) * pi$1 * x1,
      sign$1(phi) * pi$1 * sqrt$1(1 - x1 * (2 * A + x1))
    ];
  }

  vanDerGrinten2Raw.invert = function(x, y) {
    if (!x) return [0, halfPi$1 * sin$1(2 * atan$1(y / pi$1))];
    var x1 = abs$1(x / pi$1),
        A = (1 - x1 * x1 - (y /= pi$1) * y) / (2 * x1),
        A2 = A * A,
        B = sqrt$1(A2 + 1);
    return [
      sign$1(x) * pi$1 * (B - A),
      sign$1(y) * halfPi$1 * sin$1(2 * atan2$1(sqrt$1((1 - 2 * A * x1) * (A + B) - x1), sqrt$1(B + A + x1)))
    ];
  };

  function vanDerGrinten2() {
    return projection(vanDerGrinten2Raw)
        .scale(79.4183);
  }

  function vanDerGrinten3Raw(lambda, phi) {
    if (abs$1(phi) < epsilon$1) return [lambda, 0];
    var sinTheta = phi / halfPi$1,
        theta = asin$1(sinTheta);
    if (abs$1(lambda) < epsilon$1 || abs$1(abs$1(phi) - halfPi$1) < epsilon$1) return [0, pi$1 * tan$1(theta / 2)];
    var A = (pi$1 / lambda - lambda / pi$1) / 2,
        y1 = sinTheta / (1 + cos$1(theta));
    return [
      pi$1 * (sign$1(lambda) * sqrt$1(A * A + 1 - y1 * y1) - A),
      pi$1 * y1
    ];
  }

  vanDerGrinten3Raw.invert = function(x, y) {
    if (!y) return [x, 0];
    var y1 = y / pi$1,
        A = (pi$1 * pi$1 * (1 - y1 * y1) - x * x) / (2 * pi$1 * x);
    return [
      x ? pi$1 * (sign$1(x) * sqrt$1(A * A + 1) - A) : 0,
      halfPi$1 * sin$1(2 * atan$1(y1))
    ];
  };

  function vanDerGrinten3() {
    return projection(vanDerGrinten3Raw)
          .scale(79.4183);
  }

  function vanDerGrinten4Raw(lambda, phi) {
    if (!phi) return [lambda, 0];
    var phi0 = abs$1(phi);
    if (!lambda || phi0 === halfPi$1) return [0, phi];
    var B = phi0 / halfPi$1,
        B2 = B * B,
        C = (8 * B - B2 * (B2 + 2) - 5) / (2 * B2 * (B - 1)),
        C2 = C * C,
        BC = B * C,
        B_C2 = B2 + C2 + 2 * BC,
        B_3C = B + 3 * C,
        lambda0 = lambda / halfPi$1,
        lambda1 = lambda0 + 1 / lambda0,
        D = sign$1(abs$1(lambda) - halfPi$1) * sqrt$1(lambda1 * lambda1 - 4),
        D2 = D * D,
        F = B_C2 * (B2 + C2 * D2 - 1) + (1 - B2) * (B2 * (B_3C * B_3C + 4 * C2) + 12 * BC * C2 + 4 * C2 * C2),
        x1 = (D * (B_C2 + C2 - 1) + 2 * sqrt$1(F)) / (4 * B_C2 + D2);
    return [
      sign$1(lambda) * halfPi$1 * x1,
      sign$1(phi) * halfPi$1 * sqrt$1(1 + D * abs$1(x1) - x1 * x1)
    ];
  }

  vanDerGrinten4Raw.invert = function(x, y) {
    var delta;
    if (!x || !y) return [x, y];
    y /= pi$1;
    var x1 = sign$1(x) * x / halfPi$1,
        D = (x1 * x1 - 1 + 4 * y * y) / abs$1(x1),
        D2 = D * D,
        B = 2 * y,
        i = 50;
    do {
      var B2 = B * B,
          C = (8 * B - B2 * (B2 + 2) - 5) / (2 * B2 * (B - 1)),
          C_ = (3 * B - B2 * B - 10) / (2 * B2 * B),
          C2 = C * C,
          BC = B * C,
          B_C = B + C,
          B_C2 = B_C * B_C,
          B_3C = B + 3 * C,
          F = B_C2 * (B2 + C2 * D2 - 1) + (1 - B2) * (B2 * (B_3C * B_3C + 4 * C2) + C2 * (12 * BC + 4 * C2)),
          F_ = -2 * B_C * (4 * BC * C2 + (1 - 4 * B2 + 3 * B2 * B2) * (1 + C_) + C2 * (-6 + 14 * B2 - D2 + (-8 + 8 * B2 - 2 * D2) * C_) + BC * (-8 + 12 * B2 + (-10 + 10 * B2 - D2) * C_)),
          sqrtF = sqrt$1(F),
          f = D * (B_C2 + C2 - 1) + 2 * sqrtF - x1 * (4 * B_C2 + D2),
          f_ = D * (2 * C * C_ + 2 * B_C * (1 + C_)) + F_ / sqrtF - 8 * B_C * (D * (-1 + C2 + B_C2) + 2 * sqrtF) * (1 + C_) / (D2 + 4 * B_C2);
      B -= delta = f / f_;
    } while (delta > epsilon$1 && --i > 0);
    return [
      sign$1(x) * (sqrt$1(D * D + 4) + D) * pi$1 / 4,
      halfPi$1 * B
    ];
  };

  function vanDerGrinten4() {
    return projection(vanDerGrinten4Raw)
        .scale(127.16);
  }

  var A = 4 * pi$1 + 3 * sqrt$1(3),
      B = 2 * sqrt$1(2 * pi$1 * sqrt$1(3) / A);

  var wagner4Raw = mollweideBromleyRaw(B * sqrt$1(3) / pi$1, B, A / 6);

  function wagner4() {
    return projection(wagner4Raw)
        .scale(176.84);
  }

  function wagner6Raw(lambda, phi) {
    return [lambda * sqrt$1(1 - 3 * phi * phi / (pi$1 * pi$1)), phi];
  }

  wagner6Raw.invert = function(x, y) {
    return [x / sqrt$1(1 - 3 * y * y / (pi$1 * pi$1)), y];
  };

  function wagner6() {
    return projection(wagner6Raw)
        .scale(152.63);
  }

  function wagner7Raw(lambda, phi) {
    var s = 0.90631 * sin$1(phi),
        c0 = sqrt$1(1 - s * s),
        c1 = sqrt$1(2 / (1 + c0 * cos$1(lambda /= 3)));
    return [
      2.66723 * c0 * c1 * sin$1(lambda),
      1.24104 * s * c1
    ];
  }

  wagner7Raw.invert = function(x, y) {
    var t1 = x / 2.66723,
        t2 = y / 1.24104,
        p = sqrt$1(t1 * t1 + t2 * t2),
        c = 2 * asin$1(p / 2);
    return [
      3 * atan2$1(x * tan$1(c), 2.66723 * p),
      p && asin$1(y * sin$1(c) / (1.24104 * 0.90631 * p))
    ];
  };

  function wagner7() {
    return projection(wagner7Raw)
        .scale(172.632);
  }

  function wiechelRaw(lambda, phi) {
    var cosPhi = cos$1(phi),
        sinPhi = cos$1(lambda) * cosPhi,
        sin1_Phi = 1 - sinPhi,
        cosLambda = cos$1(lambda = atan2$1(sin$1(lambda) * cosPhi, -sin$1(phi))),
        sinLambda = sin$1(lambda);
    cosPhi = sqrt$1(1 - sinPhi * sinPhi);
    return [
      sinLambda * cosPhi - cosLambda * sin1_Phi,
      -cosLambda * cosPhi - sinLambda * sin1_Phi
    ];
  }

  wiechelRaw.invert = function(x, y) {
    var w = (x * x + y * y) / -2,
        k = sqrt$1(-w * (2 + w)),
        b = y * w + x * k,
        a = x * w - y * k,
        D = sqrt$1(a * a + b * b);
    return [
      atan2$1(k * b, D * (1 + w)),
      D ? -asin$1(k * a / D) : 0
    ];
  };

  function wiechel() {
    return projection(wiechelRaw)
        .rotate([0, -90, 45])
        .scale(124.75)
        .clipAngle(180 - 1e-3);
  }

  function winkel3Raw(lambda, phi) {
    var coordinates = aitoffRaw(lambda, phi);
    return [
      (coordinates[0] + lambda / halfPi$1) / 2,
      (coordinates[1] + phi) / 2
    ];
  }

  winkel3Raw.invert = function(x, y) {
    var lambda = x, phi = y, i = 25;
    do {
      var cosphi = cos$1(phi),
          sinphi = sin$1(phi),
          sin_2phi = sin$1(2 * phi),
          sin2phi = sinphi * sinphi,
          cos2phi = cosphi * cosphi,
          sinlambda = sin$1(lambda),
          coslambda_2 = cos$1(lambda / 2),
          sinlambda_2 = sin$1(lambda / 2),
          sin2lambda_2 = sinlambda_2 * sinlambda_2,
          C = 1 - cos2phi * coslambda_2 * coslambda_2,
          E = C ? acos$1(cosphi * coslambda_2) * sqrt$1(F = 1 / C) : F = 0,
          F,
          fx = 0.5 * (2 * E * cosphi * sinlambda_2 + lambda / halfPi$1) - x,
          fy = 0.5 * (E * sinphi + phi) - y,
          dxdlambda = 0.5 * F * (cos2phi * sin2lambda_2 + E * cosphi * coslambda_2 * sin2phi) + 0.5 / halfPi$1,
          dxdphi = F * (sinlambda * sin_2phi / 4 - E * sinphi * sinlambda_2),
          dydlambda = 0.125 * F * (sin_2phi * sinlambda_2 - E * sinphi * cos2phi * sinlambda),
          dydphi = 0.5 * F * (sin2phi * coslambda_2 + E * sin2lambda_2 * cosphi) + 0.5,
          denominator = dxdphi * dydlambda - dydphi * dxdlambda,
          dlambda = (fy * dxdphi - fx * dydphi) / denominator,
          dphi = (fx * dydlambda - fy * dxdlambda) / denominator;
      lambda -= dlambda, phi -= dphi;
    } while ((abs$1(dlambda) > epsilon$1 || abs$1(dphi) > epsilon$1) && --i > 0);
    return [lambda, phi];
  };

  function winkel3() {
    return projection(winkel3Raw)
        .scale(158.837);
  }



  var d3 = /*#__PURE__*/Object.freeze({
    geoAiry: airy,
    geoAiryRaw: airyRaw,
    geoAitoff: aitoff,
    geoAitoffRaw: aitoffRaw,
    geoArmadillo: armadillo,
    geoArmadilloRaw: armadilloRaw,
    geoAugust: august,
    geoAugustRaw: augustRaw,
    geoBaker: baker,
    geoBakerRaw: bakerRaw,
    geoBerghaus: berghaus,
    geoBerghausRaw: berghausRaw,
    geoBertin1953: bertin,
    geoBertin1953Raw: bertin1953Raw,
    geoBoggs: boggs,
    geoBoggsRaw: boggsRaw,
    geoBonne: bonne,
    geoBonneRaw: bonneRaw,
    geoBottomley: bottomley,
    geoBottomleyRaw: bottomleyRaw,
    geoBromley: bromley,
    geoBromleyRaw: bromleyRaw,
    geoChamberlin: chamberlin,
    geoChamberlinRaw: chamberlinRaw,
    geoChamberlinAfrica: chamberlinAfrica,
    geoCollignon: collignon,
    geoCollignonRaw: collignonRaw,
    geoCraig: craig,
    geoCraigRaw: craigRaw,
    geoCraster: craster,
    geoCrasterRaw: crasterRaw,
    geoCylindricalEqualArea: cylindricalEqualArea,
    geoCylindricalEqualAreaRaw: cylindricalEqualAreaRaw$1,
    geoCylindricalStereographic: cylindricalStereographic,
    geoCylindricalStereographicRaw: cylindricalStereographicRaw,
    geoEckert1: eckert1,
    geoEckert1Raw: eckert1Raw,
    geoEckert2: eckert2,
    geoEckert2Raw: eckert2Raw,
    geoEckert3: eckert3,
    geoEckert3Raw: eckert3Raw,
    geoEckert4: eckert4,
    geoEckert4Raw: eckert4Raw,
    geoEckert5: eckert5,
    geoEckert5Raw: eckert5Raw,
    geoEckert6: eckert6,
    geoEckert6Raw: eckert6Raw,
    geoEisenlohr: eisenlohr,
    geoEisenlohrRaw: eisenlohrRaw,
    geoFahey: fahey,
    geoFaheyRaw: faheyRaw,
    geoFoucaut: foucaut,
    geoFoucautRaw: foucautRaw,
    geoGilbert: gilbert,
    geoGingery: gingery,
    geoGingeryRaw: gingeryRaw,
    geoGinzburg4: ginzburg4,
    geoGinzburg4Raw: ginzburg4Raw,
    geoGinzburg5: ginzburg5,
    geoGinzburg5Raw: ginzburg5Raw,
    geoGinzburg6: ginzburg6,
    geoGinzburg6Raw: ginzburg6Raw,
    geoGinzburg8: ginzburg8,
    geoGinzburg8Raw: ginzburg8Raw,
    geoGinzburg9: ginzburg9,
    geoGinzburg9Raw: ginzburg9Raw,
    geoGringorten: gringorten,
    geoGringortenRaw: gringortenRaw,
    geoGuyou: guyou,
    geoGuyouRaw: guyouRaw,
    geoHammer: hammer,
    geoHammerRaw: hammerRaw,
    geoHammerRetroazimuthal: hammerRetroazimuthal,
    geoHammerRetroazimuthalRaw: hammerRetroazimuthalRaw,
    geoHealpix: healpix,
    geoHealpixRaw: healpixRaw,
    geoHill: hill,
    geoHillRaw: hillRaw,
    geoHomolosine: homolosine,
    geoHomolosineRaw: homolosineRaw,
    geoHyperelliptical: hyperelliptical,
    geoHyperellipticalRaw: hyperellipticalRaw,
    geoInterrupt: interrupt,
    geoInterruptedBoggs: boggs$1,
    geoInterruptedHomolosine: homolosine$1,
    geoInterruptedMollweide: mollweide$1,
    geoInterruptedMollweideHemispheres: mollweideHemispheres,
    geoInterruptedSinuMollweide: sinuMollweide$1,
    geoInterruptedSinusoidal: sinusoidal$1,
    geoKavrayskiy7: kavrayskiy7,
    geoKavrayskiy7Raw: kavrayskiy7Raw,
    geoLagrange: lagrange,
    geoLagrangeRaw: lagrangeRaw,
    geoLarrivee: larrivee,
    geoLarriveeRaw: larriveeRaw,
    geoLaskowski: laskowski,
    geoLaskowskiRaw: laskowskiRaw,
    geoLittrow: littrow,
    geoLittrowRaw: littrowRaw,
    geoLoximuthal: loximuthal,
    geoLoximuthalRaw: loximuthalRaw,
    geoMiller: miller,
    geoMillerRaw: millerRaw,
    geoModifiedStereographic: modifiedStereographic,
    geoModifiedStereographicRaw: modifiedStereographicRaw,
    geoModifiedStereographicAlaska: modifiedStereographicAlaska,
    geoModifiedStereographicGs48: modifiedStereographicGs48,
    geoModifiedStereographicGs50: modifiedStereographicGs50,
    geoModifiedStereographicMiller: modifiedStereographicMiller,
    geoModifiedStereographicLee: modifiedStereographicLee,
    geoMollweide: mollweide,
    geoMollweideRaw: mollweideRaw,
    geoMtFlatPolarParabolic: mtFlatPolarParabolic,
    geoMtFlatPolarParabolicRaw: mtFlatPolarParabolicRaw,
    geoMtFlatPolarQuartic: mtFlatPolarQuartic,
    geoMtFlatPolarQuarticRaw: mtFlatPolarQuarticRaw,
    geoMtFlatPolarSinusoidal: mtFlatPolarSinusoidal,
    geoMtFlatPolarSinusoidalRaw: mtFlatPolarSinusoidalRaw,
    geoNaturalEarth: naturalEarth1,
    geoNaturalEarthRaw: naturalEarth1Raw,
    geoNaturalEarth2: naturalEarth2,
    geoNaturalEarth2Raw: naturalEarth2Raw,
    geoNellHammer: nellHammer,
    geoNellHammerRaw: nellHammerRaw,
    geoPatterson: patterson,
    geoPattersonRaw: pattersonRaw,
    geoPolyconic: polyconic,
    geoPolyconicRaw: polyconicRaw,
    geoPolyhedral: polyhedral,
    geoPolyhedralButterfly: butterfly,
    geoPolyhedralCollignon: collignon$1,
    geoPolyhedralWaterman: waterman,
    geoProject: index$1,
    geoGringortenQuincuncial: gringorten$1,
    geoPeirceQuincuncial: peirce,
    geoPierceQuincuncial: peirce,
    geoQuantize: quantize,
    geoQuincuncial: quincuncial,
    geoRectangularPolyconic: rectangularPolyconic,
    geoRectangularPolyconicRaw: rectangularPolyconicRaw,
    geoRobinson: robinson,
    geoRobinsonRaw: robinsonRaw,
    geoSatellite: satellite,
    geoSatelliteRaw: satelliteRaw,
    geoSinuMollweide: sinuMollweide,
    geoSinuMollweideRaw: sinuMollweideRaw,
    geoSinusoidal: sinusoidal,
    geoSinusoidalRaw: sinusoidalRaw,
    geoStitch: stitch,
    geoTimes: times,
    geoTimesRaw: timesRaw,
    geoTwoPointAzimuthal: twoPointAzimuthal,
    geoTwoPointAzimuthalRaw: twoPointAzimuthalRaw,
    geoTwoPointAzimuthalUsa: twoPointAzimuthalUsa,
    geoTwoPointEquidistant: twoPointEquidistant,
    geoTwoPointEquidistantRaw: twoPointEquidistantRaw,
    geoTwoPointEquidistantUsa: twoPointEquidistantUsa,
    geoVanDerGrinten: vanDerGrinten,
    geoVanDerGrintenRaw: vanDerGrintenRaw,
    geoVanDerGrinten2: vanDerGrinten2,
    geoVanDerGrinten2Raw: vanDerGrinten2Raw,
    geoVanDerGrinten3: vanDerGrinten3,
    geoVanDerGrinten3Raw: vanDerGrinten3Raw,
    geoVanDerGrinten4: vanDerGrinten4,
    geoVanDerGrinten4Raw: vanDerGrinten4Raw,
    geoWagner4: wagner4,
    geoWagner4Raw: wagner4Raw,
    geoWagner6: wagner6,
    geoWagner6Raw: wagner6Raw,
    geoWagner7: wagner7,
    geoWagner7Raw: wagner7Raw,
    geoWiechel: wiechel,
    geoWiechelRaw: wiechelRaw,
    geoWinkel3: winkel3,
    geoWinkel3Raw: winkel3Raw
  });

  [
    'airy',
    'aitoff',
    'armadillo',
    'august',
    'baker',
    'berghaus',
    'bertin1953',
    'boggs',
    'bonne',
    'bottomley',
    'bromley',
    'chamberlinAfrica',
    'collignon',
    'craig',
    'craster',
    'cylindricalEqualArea',
    'cylindricalStereographic',
    'eckert1',
    'eckert2',
    'eckert3',
    'eckert4',
    'eckert5',
    'eckert6',
    'eisenlohr',
    'fahey',
    'foucaut',
    'gilbert',
    'gingery',
    'ginzburg4',
    'ginzburg5',
    'ginzburg6',
    'ginzburg8',
    'ginzburg9',
    'gringorten',
    'guyou',
    'hammer',
    'hammerRetroazimuthal',
    'healpix',
    'hill',
    'homolosine',
    'kavrayskiy7',
    'lagrange',
    'larrivee',
    'laskowski',
    'littrow',
    'loximuthal',
    'miller',
    'modifiedStereographicAlaska',
    'modifiedStereographicGs48',
    'modifiedStereographicGs50',
    'modifiedStereographicMiller',
    'modifiedStereographicLee',
    'mollweide',
    'mtFlatPolarParabolic',
    'mtFlatPolarQuartic',
    'mtFlatPolarSinusoidal',
    'naturalEarth1',
    'naturalEarth2',
    'nellHammer',
    'patterson',
    'polyconic',
    'rectangularPolyconic',
    'robinson',
    'satellite',
    'sinusoidal',
    'sinuMollweide',
    'times',
    'twoPointAzimuthalUsa',
    'twoPointEquidistantUsa',
    'vanDerGrinten',
    'vanDerGrinten2',
    'vanDerGrinten3',
    'vanDerGrinten4',
    'wagner4',
    'wagner6',
    'wagner7',
    'wiechel',
    'winkel3',
    'interruptedHomolosine',
    'interruptedSinusoidal',
    'interruptedBoggs',
    'interruptedSinuMollweide',
    'interruptedMollweide',
    'interruptedMollweideHemispheres',
    'polyhedralButterfly',
    'polyhedralCollignon',
    'polyhedralWaterman',
    'gringortenQuincuncial',
    'peirceQuincuncial'
  ].forEach(function(p) {
    vegaProjection.projection(p, d3['geo' + p[0].toUpperCase() + p.slice(1)]);
  });

})));
