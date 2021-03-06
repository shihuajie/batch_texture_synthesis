/*******************************************************************************
 * basic.h - basic patch implementation
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef BASIC_H
#define	BASIC_H

#include "../patch.h" 

#define ANGLESTEP 90
#define MIRROR 1           // define mirroring : 0 = None, 1 = MirrorX, 2 = MirrorY, 3 = MirrorXY
#define MAX_ANGLE (M_PI/2)

namespace pm {
	struct Prop {
		static const int angleCount = 2 * ANGLESTEP + 1;
		float COS[angleCount];
		float SIN[angleCount];

		Prop() {
            SIN[ANGLESTEP] = 0.0f;
			COS[ANGLESTEP] = 1.0f;
			for (int i = 1; i <= ANGLESTEP; ++i) {
				float theta = MAX_ANGLE / ANGLESTEP * i;
				float c = std::cos(theta);
                float s = std::sin(theta);
				SIN[ANGLESTEP + i] = s;
				SIN[ANGLESTEP - i] = -s;
				COS[ANGLESTEP + i] = COS[ANGLESTEP - i] = c;
			}
		}
	};
	static Prop BasicProp;
    
    template <typename T, int PatchWidth = DEFAULT_PATCH_SIZE >
	struct BasicPatch : public BasicGrid< BasicPatch<T, PatchWidth> > {
		/// the width that can be inferred statically
		static const int staticSize = PatchWidth;
		/// Location type
		typedef T Scalar;
		/// Pixel location within a texture
		typedef Point<T> PixLoc;
		typedef Point<int> Index;
		/// Patch types
		typedef BasicPatch<int, PatchWidth> OriginalPatchType;
		typedef T Coherence;
		
		// template definitions
		inline static int width(int newSize = 0) {
			if (newSize != 0) throw "Static patches cannot change their size!";
			return PatchWidth;
		}

		inline static float radius() {
			return 0.5f * width();
		}

		inline static int dimensions() {
			return 5; //modified by huajie 2015-9-23
		}

		// the real data is here
		T x, y;
        float sx, sy;
		int angleIndex; //added by huajie 2015-9-22
		// Simple constructor
		BasicPatch() : x(0), y(0),angleIndex(ANGLESTEP),sx(1),sy(1) {
		}
		BasicPatch(T y0, T x0) : x(x0), y(y0),angleIndex(ANGLESTEP),sx(1),sy(1) {
		}
		bool operator==(const BasicPatch<T, PatchWidth> &p) const {
			return x == p.x && y == p.y && angleIndex == p.angleIndex && sx == p.sx && sy == p.sy;
		}

		// points of the patch
		inline void points(PixLoc pts[4]) const {
			pts[0] = pxloc(x, y, angleIndex, sx, sy, 0, 0);
			pts[1] = pxloc(x, y, angleIndex, sx, sy, 0, width() - 1);
			pts[2] = pxloc(x, y, angleIndex, sx, sy, width() - 1, width() - 1);
			pts[3] = pxloc(x, y, angleIndex, sx, sy, width() - 1, 0);
		}

		inline int cx() const {
			return x + int(radius());
		}

		inline int cy() const {
			return y + int(radius());
		}

		template <typename Storage>
		inline void store(Storage &out, int channels) const {
			switch (channels) {
                case 5://added by huajie 2015-9-23
					out[4] = sy;
				case 4:
					out[3] = sx;
				case 3: 
					out[2] = angleIndex;
				case 2:
					out[1] = x;
					out[0] = y;
					break;
				default:
					std::cerr << "Unknown storage of " << channels;
					std::cerr << " channels!\n";
			}
		}
		
		inline void load(const float *source, int channels, int offset) { //modified by huajie 2015-9-23
			if(channels != dimensions()){
				std::cerr << "Warning: loading an incomplete source of " << channels;
				std::cerr << "channels (current=" << dimensions() << " channels)\n";
			}
			switch (channels) {
                case 5:
					sy = source[offset * 4]>0 ? 1 : -1;
				case 4:
					sx = source[offset * 3]>0 ? 1 : -1;
                case 3:
                    angleIndex = round(source[offset * 2]);
                case 2:
                    y = source[0];
                    x = source[offset];
                    break;
                default:
                    std::cerr << "Unknown source of " << channels;
                    std::cerr << " channels!\n";
			}
		}

		//addey by huajie 2015-9-22
        bool withinFrame(const Image *frame) const {
            const int maxX = frame->cols;
            const int maxY = frame->rows;
            // basics
            if (x < 0 || x > maxX - width() || y < 0 || y > maxY - width()) return false;
            // top left
            PixLoc px = pxloc(x, y, angleIndex, sx, sy, 0, 0);
            if (px.x < 0.0f || px.x >= maxX || px.y < 0.0f || px.y >= maxY) return false;
            // top right
            px = pxloc(x, y, angleIndex, sx, sy, 0, width() - 1);
            if (px.x < 0.0f || px.x >= maxX || px.y < 0.0f || px.y >= maxY) return false;
            // top left
            px = pxloc(x, y, angleIndex, sx, sy, width() - 1, 0);
            if (px.x < 0.0f || px.x >= maxX || px.y < 0.0f || px.y >= maxY) return false;
            // top left
            px = pxloc(x, y, angleIndex, sx, sy, width() - 1, width() - 1);
            if (px.x < 0.0f || px.x >= maxX || px.y < 0.0f || px.y >= maxY) return false;

            // everything is ok!
            return true;
        }

		inline static PixLoc pxloc(float x, float y, int index, float sx, float sy, int dx, int dy) {
			// translation to the center
			float px = dx;
			float py = dy;
			//float angle = (angleIndex - ANGLESTEP)*M_PI_4;
			// transformation paramaters
			float cosA = BasicProp.COS[index];
			float sinA = BasicProp.SIN[index];
			// actual transformation
			// 1 =  mirror
            px *= sx;
			py *= sy;
			// 2 = rotate
			float tx = cosA * px - sinA * py;
			float ty = sinA * px + cosA * py;
			// 3 = translate
			return PixLoc(x + tx, y + ty);
		}
		inline PixLoc transform(int py, int px) const {
			return pxloc(x, y, angleIndex, sx, sy, px, py);
		}
		inline PixLoc transform(const Index &i) const {
			return transform(i.y, i.x);
		}
		inline PixLoc operator *(const Index &i) const {
			return transform(i);
		}
		//end
	};
	template <>
	inline int BasicPatch<int, 0>::width(int newSize) {
		static int size = DEFAULT_PATCH_SIZE;
		if (newSize > 0) size = newSize;
		return size;
	}
	
	template <>
	inline int BasicPatch<float, 0>::width(int newSize) {
		return BasicPatch<int, 0>::width(newSize); // delegate
	}

	/**
	 * \brief Dynamic sized basic patch
	 */
	typedef BasicPatch<int, 0> BasicPatchX;
	typedef BasicPatch<float, 0> BasicFloatPatchX;
    #define MAX_ITERATIONS 1000
	template <typename T, int PatchWidth>
	inline void randomInit(RNG rand, const Image *parent,
			BasicPatch<T, PatchWidth> &patch) {
        typedef BasicPatch<T, PatchWidth> Patch;
		T maxX = parent->cols - Patch::width() - 1;
		T maxY = parent->rows - Patch::width() - 1;
		int it = 0;
		do {
			if (it++ > MAX_ITERATIONS) {
				std::cerr << "Spent " << it << " iterations to initialize patch (width=" ;
				std::cerr << " in [" << parent->cols << "x" << parent->rows << "])\n";
				patch.angleIndex = ANGLESTEP;
				break;
			}
			patch.x = uniform(rand, T(0), maxX);
			patch.y = uniform(rand, T(0), maxY);
			patch.angleIndex = uniform(rand, 0, 2 * ANGLESTEP);
            switch (MIRROR) {
				case 1:
					if (bernoulli(rand)) patch.sx *= -1;
					break;
				case 2:
					if (bernoulli(rand)) patch.sy *= -1;
					break;
				case 3:
					if (bernoulli(rand)) patch.sx *= -1;
					if (bernoulli(rand)) patch.sy *= -1;
					break;
				default:
					break;
			}
		} while (!patch.withinFrame(parent));
	}

	template <typename T, int PatchWidth>
	inline bool random(RNG rand, const Image *parent,
			const BasicPatch<T, PatchWidth> &oldPatch,
			BasicPatch<T, PatchWidth> &newPatch,
			int windowSize) {
		typedef BasicPatch<T, PatchWidth> Patch;
		newPatch.x = uniform(rand,
				std::max(T(0), oldPatch.x - T(windowSize)),
				std::min(T(parent->cols - Patch::width()), oldPatch.x + T(windowSize)) // Not cols - P - 1!
				);
		newPatch.y = uniform(rand,
				std::max(T(0), oldPatch.y - T(windowSize)),
				std::min(T(parent->rows - Patch::width()), oldPatch.y + T(windowSize))
				);
        newPatch.angleIndex = uniform(rand, 0, 2 * ANGLESTEP);
        switch (MIRROR) {
			case 1:
				if (bernoulli(rand)) newPatch.sx *= -1;
				break;
			case 2:
				if (bernoulli(rand)) newPatch.sy *= -1;
				break;
			case 3:
				if (bernoulli(rand)) newPatch.sx *= -1;
				if (bernoulli(rand)) newPatch.sy *= -1;
				break;
			default:
				break;
		}
		// the random patch may not be valid!
		return newPatch.withinFrame(parent);
	}
    
    template <typename T, int PatchWidth>
	inline bool aligned(RNG rand, const Image *parent,
			const BasicPatch<T, PatchWidth> &oldPatch,
            BasicPatch<T, PatchWidth> &newPatch, 
            const Point2f &g1, const Point2f &g2, float jitter) {
        typedef BasicPatch<T, PatchWidth> Patch;
		// steps and direction type
		const bool sameDir = g1.x * g1.y >= 0;
		const Point2f step = Point2f::max(g1.abs(), Point2f(1, 1));
		
		// what are the possible steps?
        int xPosSteps = std::floor((parent->rows - Patch::width() - oldPatch.x) / step.x);
		int xNegSteps = std::floor((oldPatch.x) / step.x);
		int yPosSteps = std::floor((parent->cols - Patch::width() - oldPatch.y) / step.y);
		int yNegSteps = std::floor((oldPatch.y) / step.y);
		int posSteps, negSteps;
		if(sameDir){
			// same direction in x and y
			posSteps = std::min(xPosSteps, yPosSteps);
			negSteps = std::min(xNegSteps, yNegSteps);
		} else {
			// inverse directions, we choose arbitrarily the association
			posSteps = std::min(xPosSteps, yNegSteps);
			negSteps = std::min(xNegSteps, yPosSteps);
		}
		
		// can we go anywhere?
		if(posSteps + negSteps < 1){
			// no transformation
			newPatch.x = oldPatch.x;
			newPatch.y = oldPatch.y;
			// we failed
			return false;
		}
		
		// get our step
		int stepID = uniform(rand, 0, posSteps + negSteps);
		// direction
		int dir;
		if(stepID > posSteps){
			stepID = 1 + stepID - posSteps;
			dir = -1;
		} else if(posSteps) {
			dir = 1;
		}
		
		// compute the shift
		Point2f shift;
		if(sameDir){
			shift = (g1.x > 0 || g1.y > 0 ? g1 : -g1) * (dir * stepID);
		} else if(g1.x > 0) { // => g1.y < 0
			shift = g1 * (dir * stepID);
		} else { // => g1.y > 0
			shift = -g1 * (dir * stepID); 
			// reversed because of the definition of posSteps and negSteps
			// i.e. taken for g1.x positive
		}
		
		// is it the end or do we have to go for g2?
		if(g2.isOrigin()) {
			if(jitter == 0.0f) {
				newPatch.x = oldPatch.x + roundOrNot<T>(shift.x);
				newPatch.y = oldPatch.y + roundOrNot<T>(shift.y);
			} else {
				// add the jitter from a gaussian
				T dx = gaussian(rand, jitter); // /!\ dx and dy may be correlated, but that's fine
				T dy = gaussian(rand, jitter);
				Patch jittered(
					roundOrNot<T>(oldPatch.x + shift.x + dx),
					roundOrNot<T>(oldPatch.x + shift.y + dy)
				);
				if(isWithin(parent, jittered)){
					newPatch = jittered;
				}
			}
		} else {
			// we go for g2, recursively of course ...
			Patch tmpPatch(oldPatch.x + shift.x, oldPatch.y + shift.y);
			aligned(rand, parent, tmpPatch, newPatch, g2, Point2f(), jitter);
		}
		return isWithin(parent, newPatch);
    }

	template <typename T, int PatchWidth>
	inline void deltaPatch(
			const BasicPatch<T, PatchWidth> &patch,
			BasicPatch<T, PatchWidth> &delta,
			int dy, int dx) {
 				typedef BasicPatch<T, PatchWidth> Patch;
				typedef typename Patch::PixLoc PixLoc;
				delta = patch;
				PixLoc dp = Patch::pxloc(patch.x, patch.y, patch.angleIndex, patch.sx, patch.sy, dx, dy); //modified by huajie 2015-9-23
				delta.x = dp.x;
				delta.y = dp.y;
	}

	template <typename T, int PatchWidth>
	inline bool isWithin(const Image *parent, const BasicPatch<T, PatchWidth> &patch) {
		return patch.withinFrame(parent);
	}
	template <int PatchWidth>
	inline typename BasicPatch<int, PatchWidth>::Coherence coherence(
			const BasicPatch<int, PatchWidth> &p1, 
			const BasicPatch<int, PatchWidth> &p2, 
			int dy, int dx){
		return  p1.y + dy == p2.y && p1.x + dx == p2.x ? 1 : 0;
	}
	template <int PatchWidth>
	inline typename BasicPatch<float, PatchWidth>::Coherence coherence(
			const BasicPatch<float, PatchWidth> &p1, 
			const BasicPatch<float, PatchWidth> &p2, 
			int dy, int dx){
		return p1.y + dy == p2.y && p1.x + dx == p2.x ? 1.0f : 0.0f;
		// return  std::max(0.0f, 1.0f - std::abs(p1.y + dy - p2.y) - std::abs(p1.x + dx - p2.x));
	}
}

#endif	/* BASIC_H */