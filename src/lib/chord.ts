/**
 * Chord representation using dual bitmask + Center of Gravity (CoG)
 *
 * Bitmask: 12-bit integer where bit N = semitone N is present (C=0, C#=1, ..., B=11)
 * CoG: 3D coordinate based on diminished tetrad axes
 *
 * This gives us O(1) equality/hashing AND geometric meaning.
 */

// ============================================================================
// CONSTANTS
// ============================================================================

/** Note names in chromatic order */
export const CHROMATIC = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B'] as const;
export type NoteName = typeof CHROMATIC[number];

/** Diminished tetrad axes - each axis contains 4 notes a minor 3rd apart */
export const AXIS = {
	X: ['C', 'Eb', 'Gb', 'A'] as const,   // Dim7 built on C
	Y: ['D', 'F', 'Ab', 'B'] as const,    // Dim7 built on D
	Z: ['E', 'G', 'Bb', 'Db'] as const,   // Dim7 built on E
} as const;

/** Semitone index for each note (C=0) */
export const SEMI: Record<NoteName, number> = {
	C: 0, Db: 1, D: 2, Eb: 3, E: 4, F: 5,
	Gb: 6, G: 7, Ab: 8, A: 9, Bb: 10, B: 11
};

/** Reverse lookup: semitone → note name */
export const NOTE_FROM_SEMI: NoteName[] = CHROMATIC.slice();

/** Which axis each note belongs to (0=X, 1=Y, 2=Z) */
export const NOTE_AXIS: Record<NoteName, 0 | 1 | 2> = {
	C: 0, Eb: 0, Gb: 0, A: 0,   // X axis
	D: 1, F: 1, Ab: 1, B: 1,    // Y axis
	E: 2, G: 2, Bb: 2, Db: 2,   // Z axis
};

// ============================================================================
// BITMASK OPERATIONS
// ============================================================================

/** Convert note names to a 12-bit bitmask */
export function notesToBits(notes: NoteName[]): number {
	let bits = 0;
	for (const note of notes) {
		bits |= (1 << SEMI[note]);
	}
	return bits;
}

/** Convert semitone indices to a 12-bit bitmask */
export function semitonesToBits(semitones: number[]): number {
	let bits = 0;
	for (const semi of semitones) {
		bits |= (1 << (semi % 12));
	}
	return bits;
}

/** Convert a bitmask back to note names */
export function bitsToNotes(bits: number): NoteName[] {
	const notes: NoteName[] = [];
	for (let i = 0; i < 12; i++) {
		if (bits & (1 << i)) {
			notes.push(NOTE_FROM_SEMI[i]);
		}
	}
	return notes;
}

/** Convert a bitmask to semitone indices */
export function bitsToSemitones(bits: number): number[] {
	const semis: number[] = [];
	for (let i = 0; i < 12; i++) {
		if (bits & (1 << i)) {
			semis.push(i);
		}
	}
	return semis;
}

/** Count the number of notes in a bitmask (popcount) */
export function bitCount(bits: number): number {
	// Brian Kernighan's algorithm
	let count = 0;
	while (bits) {
		bits &= bits - 1;
		count++;
	}
	return count;
}

/**
 * Transpose a chord by N semitones (rotate bits with wrapping)
 * Positive = up, Negative = down
 */
export function transpose(bits: number, semitones: number): number {
	// Normalize semitones to 0-11 range
	const shift = ((semitones % 12) + 12) % 12;
	if (shift === 0) return bits;

	// Rotate left by shift, wrapping at 12 bits
	const shifted = bits << shift;
	return ((shifted | (shifted >> 12)) & 0xFFF);
}

/** Get the complement of a chord (all notes NOT in the chord) */
export function complement(bits: number): number {
	return (~bits) & 0xFFF;
}

/** Union of two chords (all notes in either) */
export function union(a: number, b: number): number {
	return a | b;
}

/** Intersection of two chords (notes in both) */
export function intersection(a: number, b: number): number {
	return a & b;
}

/** Symmetric difference (notes in one but not both) */
export function symmetricDiff(a: number, b: number): number {
	return a ^ b;
}

/** Check if chord A is a subset of chord B */
export function isSubset(a: number, b: number): boolean {
	return (a & b) === a;
}

/** Check if chord A is a superset of chord B */
export function isSuperset(a: number, b: number): boolean {
	return (a & b) === b;
}

// ============================================================================
// CENTER OF GRAVITY (CoG) OPERATIONS
// ============================================================================

export type Coord3D = [number, number, number];

/** Calculate CoG coordinate from a bitmask */
export function bitsToCoG(bits: number): Coord3D {
	let x = 0, y = 0, z = 0;
	for (let i = 0; i < 12; i++) {
		if (bits & (1 << i)) {
			const axis = NOTE_AXIS[NOTE_FROM_SEMI[i]];
			if (axis === 0) x++;
			else if (axis === 1) y++;
			else z++;
		}
	}
	return [x, y, z];
}

/** Calculate CoG coordinate from note names */
export function notesToCoG(notes: NoteName[]): Coord3D {
	let x = 0, y = 0, z = 0;
	for (const note of notes) {
		const axis = NOTE_AXIS[note];
		if (axis === 0) x++;
		else if (axis === 1) y++;
		else z++;
	}
	return [x, y, z];
}

/** Get the complement coordinate: n-note at (x,y,z) → (12-n)-note at (4-x, 4-y, 4-z) */
export function complementCoord(coord: Coord3D): Coord3D {
	return [4 - coord[0], 4 - coord[1], 4 - coord[2]];
}

/** Euclidean distance between two CoG points */
export function cogDistance(a: Coord3D, b: Coord3D): number {
	const dx = a[0] - b[0];
	const dy = a[1] - b[1];
	const dz = a[2] - b[2];
	return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/** Manhattan distance between two CoG points */
export function cogManhattan(a: Coord3D, b: Coord3D): number {
	return Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]) + Math.abs(a[2] - b[2]);
}

/** Convert Coord3D to a string key for Map/Set usage */
export function coordKey(coord: Coord3D): string {
	return `${coord[0]},${coord[1]},${coord[2]}`;
}

// ============================================================================
// INTERVAL OPERATIONS
// ============================================================================

/** Get all intervals present in a chord (as semitone distances) */
export function getIntervals(bits: number): number[] {
	const semis = bitsToSemitones(bits);
	const intervals = new Set<number>();

	for (let i = 0; i < semis.length; i++) {
		for (let j = i + 1; j < semis.length; j++) {
			const interval = semis[j] - semis[i];
			intervals.add(interval);
		}
	}

	return Array.from(intervals).sort((a, b) => a - b);
}

/** Get the interval vector (how many of each interval class 1-6) */
export function intervalVector(bits: number): [number, number, number, number, number, number] {
	const semis = bitsToSemitones(bits);
	const vector: [number, number, number, number, number, number] = [0, 0, 0, 0, 0, 0];

	for (let i = 0; i < semis.length; i++) {
		for (let j = i + 1; j < semis.length; j++) {
			let interval = semis[j] - semis[i];
			// Fold to interval class (1-6)
			if (interval > 6) interval = 12 - interval;
			if (interval >= 1 && interval <= 6) {
				vector[interval - 1]++;
			}
		}
	}

	return vector;
}

// ============================================================================
// CHORD CLASS - THE MAIN DATA STRUCTURE
// ============================================================================

export class Chord {
	/** 12-bit bitmask representation */
	readonly bits: number;

	/** Center of Gravity in the diminished axis coordinate system */
	readonly cog: Coord3D;

	/** Cached interval vector */
	private _intervalVector?: [number, number, number, number, number, number];

	constructor(bits: number) {
		this.bits = bits & 0xFFF; // Ensure 12-bit
		this.cog = bitsToCoG(this.bits);
	}

	// --- Static constructors ---

	static fromNotes(notes: NoteName[]): Chord {
		return new Chord(notesToBits(notes));
	}

	static fromSemitones(semitones: number[]): Chord {
		return new Chord(semitonesToBits(semitones));
	}

	static fromIntervals(root: number, intervals: number[]): Chord {
		const semitones = [root, ...intervals.map(i => (root + i) % 12)];
		return new Chord(semitonesToBits(semitones));
	}

	// --- Getters ---

	get notes(): NoteName[] {
		return bitsToNotes(this.bits);
	}

	get semitones(): number[] {
		return bitsToSemitones(this.bits);
	}

	get size(): number {
		return bitCount(this.bits);
	}

	get intervalVector(): [number, number, number, number, number, number] {
		if (!this._intervalVector) {
			this._intervalVector = intervalVector(this.bits);
		}
		return this._intervalVector;
	}

	// --- Transformations (return new Chord) ---

	transpose(semitones: number): Chord {
		return new Chord(transpose(this.bits, semitones));
	}

	complement(): Chord {
		return new Chord(complement(this.bits));
	}

	union(other: Chord): Chord {
		return new Chord(union(this.bits, other.bits));
	}

	intersection(other: Chord): Chord {
		return new Chord(intersection(this.bits, other.bits));
	}

	// --- Comparisons ---

	equals(other: Chord): boolean {
		return this.bits === other.bits;
	}

	isSubsetOf(other: Chord): boolean {
		return isSubset(this.bits, other.bits);
	}

	isSupersetOf(other: Chord): boolean {
		return isSuperset(this.bits, other.bits);
	}

	/** Distance in CoG space */
	distanceTo(other: Chord): number {
		return cogDistance(this.cog, other.cog);
	}

	/** Check if has specific interval (in semitones) */
	hasInterval(interval: number): boolean {
		const semis = this.semitones;
		for (let i = 0; i < semis.length; i++) {
			for (let j = i + 1; j < semis.length; j++) {
				if (semis[j] - semis[i] === interval) return true;
				if (semis[j] - semis[i] === 12 - interval) return true;
			}
		}
		return false;
	}

	/** Check for perfect fifth (7 semitones) */
	hasPerfectFifth(): boolean {
		return this.hasInterval(7);
	}

	/** Check for perfect fourth (5 semitones) */
	hasPerfectFourth(): boolean {
		return this.hasInterval(5);
	}

	/** Check for tritone (6 semitones) */
	hasTritone(): boolean {
		return this.hasInterval(6);
	}

	// --- String representation ---

	toString(): string {
		return this.notes.join(' ');
	}

	toDebug(): string {
		return `Chord(${this.bits.toString(2).padStart(12, '0')}) = [${this.notes.join(', ')}] @ (${this.cog.join(', ')})`;
	}
}

// ============================================================================
// PRECOMPUTED LOOKUP TABLE
// ============================================================================

export interface ChordEntry {
	chord: Chord;
	name?: string;
}

/**
 * Precomputed lookup table for all 4096 possible pitch class sets.
 * Since we're using bitmasks, we can index directly by the integer value.
 */
export class ChordLookup {
	private table: ChordEntry[] = new Array(4096);

	constructor() {
		// Precompute all 4096 chords
		for (let bits = 0; bits < 4096; bits++) {
			this.table[bits] = {
				chord: new Chord(bits),
				name: this.identifyChordType(bits)
			};
		}
	}

	/** Get chord entry by bitmask (O(1)) */
	get(bits: number): ChordEntry {
		return this.table[bits & 0xFFF];
	}

	/** Get chord entry by notes */
	getByNotes(notes: NoteName[]): ChordEntry {
		return this.get(notesToBits(notes));
	}

	/** Get all chords at a specific CoG coordinate */
	getAtCoord(coord: Coord3D): ChordEntry[] {
		const key = coordKey(coord);
		return this.table.filter(entry => coordKey(entry.chord.cog) === key);
	}

	/** Get all chords of a specific size */
	getBySize(size: number): ChordEntry[] {
		return this.table.filter(entry => entry.chord.size === size);
	}

	/** Find chords by name pattern */
	findByName(pattern: string | RegExp): ChordEntry[] {
		const regex = typeof pattern === 'string' ? new RegExp(pattern, 'i') : pattern;
		return this.table.filter(entry => entry.name && regex.test(entry.name));
	}

	/** Group all chords by their CoG coordinate */
	groupByCoord(): Map<string, ChordEntry[]> {
		const groups = new Map<string, ChordEntry[]>();
		for (const entry of this.table) {
			const key = coordKey(entry.chord.cog);
			if (!groups.has(key)) {
				groups.set(key, []);
			}
			groups.get(key)!.push(entry);
		}
		return groups;
	}

	/** Identify common chord types */
	private identifyChordType(bits: number): string | undefined {
		const semis = bitsToSemitones(bits);
		if (semis.length < 2) return semis.length === 1 ? 'Note' : 'Empty';

		// Get intervals from lowest note
		const fromRoot = semis.slice(1).map(s => s - semis[0]);

		// Check all rotations for standard chord patterns
		for (let rotation = 0; rotation < semis.length; rotation++) {
			const rotated = rotation === 0
				? fromRoot
				: this.rotateIntervals(semis, rotation);

			const name = this.matchPattern(rotated);
			if (name) {
				return rotation === 0 ? name : `${name} inv`;
			}
		}

		return `${semis.length}-note`;
	}

	private rotateIntervals(semis: number[], rotation: number): number[] {
		const rotated = [...semis.slice(rotation), ...semis.slice(0, rotation).map(s => s + 12)];
		return rotated.slice(1).map(s => s - rotated[0]);
	}

	private matchPattern(intervals: number[]): string | undefined {
		const patterns: Record<string, number[]> = {
			// Triads
			'Major': [4, 7],
			'Minor': [3, 7],
			'Dim': [3, 6],
			'Aug': [4, 8],
			'Sus4': [5, 7],
			'Sus2': [2, 7],
			// 7th chords
			'Maj7': [4, 7, 11],
			'7': [4, 7, 10],
			'Min7': [3, 7, 10],
			'MinMaj7': [3, 7, 11],
			'Dim7': [3, 6, 9],
			'ø7': [3, 6, 10],
			'Aug7': [4, 8, 10],
			// Intervals
			'Min2': [1],
			'Maj2': [2],
			'Min3': [3],
			'Maj3': [4],
			'P4': [5],
			'Tritone': [6],
			'P5': [7],
			'Min6': [8],
			'Maj6': [9],
			'Min7Int': [10],
			'Maj7Int': [11],
		};

		for (const [name, pattern] of Object.entries(patterns)) {
			if (intervals.length === pattern.length &&
				intervals.every((v, i) => v === pattern[i])) {
				return name;
			}
		}
		return undefined;
	}
}

// ============================================================================
// CONVENIENCE EXPORTS
// ============================================================================

/** Singleton lookup table (lazy initialized) */
let _lookup: ChordLookup | null = null;
export function getLookup(): ChordLookup {
	if (!_lookup) {
		_lookup = new ChordLookup();
	}
	return _lookup;
}

// Common chord shortcuts
export const COMMON_CHORDS = {
	// Major triads
	C: Chord.fromNotes(['C', 'E', 'G']),
	Db: Chord.fromNotes(['Db', 'F', 'Ab']),
	D: Chord.fromNotes(['D', 'Gb', 'A']),
	Eb: Chord.fromNotes(['Eb', 'G', 'Bb']),
	E: Chord.fromNotes(['E', 'Ab', 'B']),
	F: Chord.fromNotes(['F', 'A', 'C']),
	Gb: Chord.fromNotes(['Gb', 'Bb', 'Db']),
	G: Chord.fromNotes(['G', 'B', 'D']),
	Ab: Chord.fromNotes(['Ab', 'C', 'Eb']),
	A: Chord.fromNotes(['A', 'Db', 'E']),
	Bb: Chord.fromNotes(['Bb', 'D', 'F']),
	B: Chord.fromNotes(['B', 'Eb', 'Gb']),

	// Minor triads
	Cm: Chord.fromNotes(['C', 'Eb', 'G']),
	Dm: Chord.fromNotes(['D', 'F', 'A']),
	Em: Chord.fromNotes(['E', 'G', 'B']),
	Fm: Chord.fromNotes(['F', 'Ab', 'C']),
	Gm: Chord.fromNotes(['G', 'Bb', 'D']),
	Am: Chord.fromNotes(['A', 'C', 'E']),
	Bm: Chord.fromNotes(['B', 'D', 'Gb']),

	// Diminished 7ths (the three axes!)
	DimC: Chord.fromNotes(['C', 'Eb', 'Gb', 'A']),
	DimD: Chord.fromNotes(['D', 'F', 'Ab', 'B']),
	DimE: Chord.fromNotes(['E', 'G', 'Bb', 'Db']),
} as const;
