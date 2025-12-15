<script lang="ts">
	import { onMount } from 'svelte';
	import * as THREE from 'three';
	import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
	import {
		Chord,
		getLookup,
		notesToBits,
		bitsToNotes,
		transpose,
		complement,
		type NoteName,
		SEMI
	} from '$lib/chord';

	// All 12 notes assigned to their diminished tetrad axis
	const AXIS_X = ['C', 'Eb', 'Gb', 'A']; // Dim C
	const AXIS_Y = ['D', 'F', 'Ab', 'B']; // Dim D
	const AXIS_Z = ['E', 'G', 'Bb', 'Db']; // Dim E
	const ALL_NOTES = [...AXIS_X, ...AXIS_Y, ...AXIS_Z];

	const AXIS_X_SET = new Set(AXIS_X);
	const AXIS_Y_SET = new Set(AXIS_Y);
	const AXIS_Z_SET = new Set(AXIS_Z);

	// Generate all combinations of n notes from the 12
	function combinations<T>(arr: T[], n: number): T[][] {
		if (n === 0) return [[]];
		if (arr.length < n) return [];
		const [first, ...rest] = arr;
		const withFirst = combinations(rest, n - 1).map((c) => [first, ...c]);
		const withoutFirst = combinations(rest, n);
		return [...withFirst, ...withoutFirst];
	}

	// Calculate coordinate for a set of notes
	function notesToCoord(notes: string[]): [number, number, number] {
		let x = 0,
			y = 0,
			z = 0;
		for (const note of notes) {
			if (AXIS_X_SET.has(note)) x++;
			else if (AXIS_Y_SET.has(note)) y++;
			else if (AXIS_Z_SET.has(note)) z++;
		}
		return [x, y, z];
	}

	// Note to semitone mapping (C = 0)
	const NOTE_TO_SEMI: Record<string, number> = {
		C: 0,
		Db: 1,
		D: 2,
		Eb: 3,
		E: 4,
		F: 5,
		Gb: 6,
		G: 7,
		Ab: 8,
		A: 9,
		Bb: 10,
		B: 11
	};

	// ========== BITMASK DEMO STATE ==========
	let showBitmaskDemo = $state(false);
	let demoChord = $state(Chord.fromNotes(['C', 'E', 'G']));
	let demoTranspose = $state(0);

	// Get all intervals present in a chord (as semitone counts)
	function getIntervals(notes: string[]): number[] {
		const semis = notes.map((n) => NOTE_TO_SEMI[n]).sort((a, b) => a - b);
		const intervals: Set<number> = new Set();

		for (let i = 0; i < semis.length; i++) {
			for (let j = i + 1; j < semis.length; j++) {
				const interval = (semis[j] - semis[i]) % 12;
				intervals.add(interval);
				// Also add the complement (for inversional equivalence)
				intervals.add((12 - interval) % 12);
			}
		}
		return Array.from(intervals).sort((a, b) => a - b);
	}

	// Check for presence of Perfect 4th (5) or Perfect 5th (7)
	type FourthFifthStatus = 'P4' | 'P5' | 'both' | 'none';
	function hasFourthOrFifth(notes: string[]): FourthFifthStatus {
		const intervals = getIntervals(notes);
		const hasP4 = intervals.includes(5);
		const hasP5 = intervals.includes(7);
		if (hasP4 && hasP5) return 'both';
		if (hasP4) return 'P4';
		if (hasP5) return 'P5';
		return 'none';
	}

	// Check for presence of minor 2nd (1) or Major 7th (11)
	type SemitoneMaj7Status = 'm2' | 'M7' | 'both' | 'none';
	function hasSemitoneMaj7(notes: string[]): SemitoneMaj7Status {
		const intervals = getIntervals(notes);
		const hasM2 = intervals.includes(1);
		const hasMaj7 = intervals.includes(11);
		if (hasM2 && hasMaj7) return 'both';
		if (hasM2) return 'm2';
		if (hasMaj7) return 'M7';
		return 'none';
	}

	// Map note to its diminished axis (0=X/DimC, 1=Y/DimD, 2=Z/DimE)
	function noteToAxis(note: string): number {
		if (AXIS_X_SET.has(note)) return 0;
		if (AXIS_Y_SET.has(note)) return 1;
		return 2;
	}

	// Calculate interval pattern: how intervals distribute across diminished axes
	// Returns counts for each axis transition type
	// ↑ = within same axis (0), ← = X→Y or Y→Z or Z→X (+1 axis), ↓ = complement, → = X→Z etc (-1 axis)
	type IntervalPattern = {
		same: number; // intervals within same diminished axis (↑)
		forward: number; // intervals moving +1 axis (←)
		backward: number; // intervals moving -1 axis (→)
		arrows: string; // visual representation
	};

	function getIntervalPattern(notes: string[]): IntervalPattern {
		const semis = notes.map((n) => NOTE_TO_SEMI[n]).sort((a, b) => a - b);
		let same = 0,
			forward = 0,
			backward = 0;

		// For each adjacent semitone pair in sorted order
		for (let i = 0; i < semis.length; i++) {
			const note1 = Object.entries(NOTE_TO_SEMI).find(([, v]) => v === semis[i])![0];
			const note2Idx = (i + 1) % semis.length;
			const note2 = Object.entries(NOTE_TO_SEMI).find(([, v]) => v === semis[note2Idx])![0];

			const axis1 = noteToAxis(note1);
			const axis2 = noteToAxis(note2);

			const diff = (axis2 - axis1 + 3) % 3;
			if (diff === 0) same++;
			else if (diff === 1) forward++;
			else backward++;
		}

		// Create arrow string (sorted for consistency)
		const arrows =
			'↑'.repeat(same) + '←'.repeat(forward) + '→'.repeat(backward);

		return { same, forward, backward, arrows };
	}

	// Identify chord type from notes
	function identifyChordType(notes: string[]): string {
		const semis = notes.map((n) => NOTE_TO_SEMI[n]).sort((a, b) => a - b);
		if (semis.length < 2) return 'Note';

		// Get intervals from root
		const fromRoot = semis.slice(1).map((s) => s - semis[0]);

		// Common patterns (intervals from root)
		const patterns: Record<string, number[]> = {
			// Triads
			Major: [4, 7],
			Minor: [3, 7],
			Dim: [3, 6],
			Aug: [4, 8],
			Sus4: [5, 7],
			Sus2: [2, 7],
			// 7th chords
			Maj7: [4, 7, 11],
			'7': [4, 7, 10],
			Min7: [3, 7, 10],
			MinMaj7: [3, 7, 11],
			Dim7: [3, 6, 9],
			'ø7': [3, 6, 10],
			Aug7: [4, 8, 10],
			// Intervals
			Min2: [1],
			Maj2: [2],
			Min3: [3],
			Maj3: [4],
			P4: [5],
			Tritone: [6],
			P5: [7],
			Min6: [8],
			Maj6: [9],
			Min7Int: [10],
			Maj7Int: [11]
		};

		// Check for matches (need to check all rotations for chords)
		for (const [name, pattern] of Object.entries(patterns)) {
			if (fromRoot.length === pattern.length) {
				if (fromRoot.every((v, i) => v === pattern[i])) {
					return name;
				}
			}
		}

		// Check rotations for triads/7ths
		if (semis.length === 3 || semis.length === 4) {
			for (let rotation = 1; rotation < semis.length; rotation++) {
				const rotated = [...semis.slice(rotation), ...semis.slice(0, rotation).map((s) => s + 12)];
				const rotatedFromRoot = rotated.slice(1).map((s) => s - rotated[0]);

				for (const [name, pattern] of Object.entries(patterns)) {
					if (rotatedFromRoot.length === pattern.length) {
						if (rotatedFromRoot.every((v, i) => v === pattern[i])) {
							return name + ' inv';
						}
					}
				}
			}
		}

		return `${semis.length}-note`;
	}

	// Generate all possible chords grouped by size
	type ChordData = {
		notes: string[];
		coord: [number, number, number];
		chordType?: string;
		fourthFifth: FourthFifthStatus;
		semitoneMaj7: SemitoneMaj7Status;
		intervalPattern: IntervalPattern;
	};
	type CoordPoint = {
		coord: [number, number, number];
		chords: ChordData[];
		count: number;
	};

	function generateAllChords(): Map<number, Map<string, CoordPoint>> {
		const result = new Map<number, Map<string, CoordPoint>>();

		for (let size = 1; size <= 12; size++) {
			const coordMap = new Map<string, CoordPoint>();
			const combos = combinations(ALL_NOTES, size);

			for (const notes of combos) {
				const coord = notesToCoord(notes);
				const key = coord.join(',');

				if (!coordMap.has(key)) {
					coordMap.set(key, { coord, chords: [], count: 0 });
				}
				const point = coordMap.get(key)!;
				const chordType = identifyChordType(notes);
				const fourthFifth = hasFourthOrFifth(notes);
				const semitoneMaj7 = hasSemitoneMaj7(notes);
				const intervalPattern = getIntervalPattern(notes);
				point.chords.push({ notes, coord, chordType, fourthFifth, semitoneMaj7, intervalPattern });
				point.count++;
			}

			result.set(size, coordMap);
		}

		return result;
	}

	const allChordsBySize = generateAllChords();

	// Check if two chords differ by exactly one semitone (one note moved up or down by 1)
	function areSemitoneAdjacent(notes1: string[], notes2: string[]): boolean {
		if (notes1.length !== notes2.length) return false;

		const set1 = new Set(notes1);
		const set2 = new Set(notes2);

		// Find notes that differ
		const onlyIn1 = notes1.filter((n) => !set2.has(n));
		const onlyIn2 = notes2.filter((n) => !set1.has(n));

		// Must have exactly one note different
		if (onlyIn1.length !== 1 || onlyIn2.length !== 1) return false;

		// Check if they're a semitone apart
		const semi1 = NOTE_TO_SEMI[onlyIn1[0]];
		const semi2 = NOTE_TO_SEMI[onlyIn2[0]];
		const diff = Math.abs(semi1 - semi2);

		return diff === 1 || diff === 11; // 11 = wrapping around (B to C)
	}

	// Get all semitone-adjacent pairs of coordinates for a given size
	function getSemitoneEdges(size: number): Array<[[number, number, number], [number, number, number]]> {
		const coordMap = allChordsBySize.get(size)!;
		const edges: Array<[[number, number, number], [number, number, number]]> = [];
		const seen = new Set<string>();

		const points = Array.from(coordMap.values());

		for (let i = 0; i < points.length; i++) {
			for (let j = i + 1; j < points.length; j++) {
				const p1 = points[i];
				const p2 = points[j];

				// Check if any chord at p1 is semitone-adjacent to any chord at p2
				let hasEdge = false;
				for (const c1 of p1.chords) {
					if (hasEdge) break;
					for (const c2 of p2.chords) {
						if (areSemitoneAdjacent(c1.notes, c2.notes)) {
							hasEdge = true;
							break;
						}
					}
				}

				if (hasEdge) {
					const key = [p1.coord.join(','), p2.coord.join(',')].sort().join('|');
					if (!seen.has(key)) {
						seen.add(key);
						edges.push([p1.coord, p2.coord]);
					}
				}
			}
		}

		return edges;
	}

	// Colors for each chord size (rainbow gradient)
	const SIZE_COLORS: Record<number, number> = {
		1: 0xff6b6b, // Single notes - red
		2: 0xffa502, // Intervals - orange
		3: 0xffd93d, // Triads - yellow
		4: 0x6bcb77, // 7th chords - green
		5: 0x4d96ff, // 9th chords - blue
		6: 0x9b59b6, // 11th chords - purple
		7: 0xe91e63, // 13th chords - pink
		8: 0x00bcd4, // 8 notes - cyan
		9: 0x8bc34a, // 9 notes - lime
		10: 0xff5722, // 10 notes - deep orange
		11: 0x607d8b, // 11 notes - blue grey
		12: 0xffffff // All 12 - white
	};

	let container: HTMLDivElement;
	let visibleSizes: Set<number> = $state(new Set([3])); // Start with triads visible
	let showEdges: boolean = $state(true); // Show semitone adjacency edges
	let showComplements: boolean = $state(false); // Show complement connections
	let scene: THREE.Scene;
	let sizeGroups: Map<number, THREE.Group> = new Map();
	let edgeGroups: Map<number, THREE.Group> = new Map();
	let complementGroup: THREE.Group | null = null;
	let selectedPoint: CoordPoint | null = $state(null);
	let raycaster: THREE.Raycaster;
	let mouse: THREE.Vector2;
	let camera: THREE.PerspectiveCamera;
	let pointToData: Map<THREE.Object3D, { size: number; point: CoordPoint }> = new Map();

	function createTextSprite(text: string, color: number): THREE.Sprite {
		const canvas = document.createElement('canvas');
		const ctx = canvas.getContext('2d')!;
		canvas.width = 256;
		canvas.height = 64;

		ctx.fillStyle = '#' + color.toString(16).padStart(6, '0');
		ctx.font = 'bold 24px Arial';
		ctx.textAlign = 'center';
		ctx.fillText(text, 128, 40);

		const texture = new THREE.CanvasTexture(canvas);
		const material = new THREE.SpriteMaterial({ map: texture, transparent: true });
		const sprite = new THREE.Sprite(material);
		sprite.scale.set(1, 0.25, 1);
		return sprite;
	}

	function addSizeGroupToScene(size: number) {
		const coordMap = allChordsBySize.get(size)!;
		const group = new THREE.Group();
		const color = SIZE_COLORS[size];

		for (const [, point] of coordMap) {
			const [x, y, z] = point.coord;

			// Size sphere based on count (more chords at this point = bigger sphere)
			const radius = 0.08 + Math.log(point.count + 1) * 0.04;
			const geometry = new THREE.SphereGeometry(radius, 16, 16);
			const material = new THREE.MeshBasicMaterial({
				color,
				transparent: true,
				opacity: 0.85
			});
			const sphere = new THREE.Mesh(geometry, material);
			// Remap: X→X (right), Z→Y (up), Y→Z (depth)
			sphere.position.set(x, z, y);
			group.add(sphere);

			// Store reference for raycasting
			pointToData.set(sphere, { size, point });
		}

		scene.add(group);
		sizeGroups.set(size, group);

		// Add edges for semitone adjacency
		if (showEdges) {
			addEdgesForSize(size);
		}
	}

	function addEdgesForSize(size: number) {
		const edges = getSemitoneEdges(size);
		const edgeGroup = new THREE.Group();
		const color = SIZE_COLORS[size];

		for (const [coord1, coord2] of edges) {
			const [x1, y1, z1] = coord1;
			const [x2, y2, z2] = coord2;

			const lineGeometry = new THREE.BufferGeometry().setFromPoints([
				new THREE.Vector3(x1, z1, y1), // Remap: X→X, Z→Y, Y→Z
				new THREE.Vector3(x2, z2, y2)
			]);
			const lineMaterial = new THREE.LineBasicMaterial({
				color,
				transparent: true,
				opacity: 0.4
			});
			const line = new THREE.Line(lineGeometry, lineMaterial);
			edgeGroup.add(line);
		}

		scene.add(edgeGroup);
		edgeGroups.set(size, edgeGroup);
	}

	function removeEdgesForSize(size: number) {
		const edgeGroup = edgeGroups.get(size);
		if (edgeGroup) {
			scene.remove(edgeGroup);
			edgeGroups.delete(size);
		}
	}

	function removeSizeGroupFromScene(size: number) {
		const group = sizeGroups.get(size);
		if (group) {
			// Remove point references
			group.traverse((obj: THREE.Object3D) => {
				pointToData.delete(obj);
			});
			scene.remove(group);
			sizeGroups.delete(size);
		}
		// Also remove edges
		removeEdgesForSize(size);
	}

	function toggleSize(size: number) {
		if (visibleSizes.has(size)) {
			visibleSizes.delete(size);
			removeSizeGroupFromScene(size);
		} else {
			visibleSizes.add(size);
			addSizeGroupToScene(size);
		}
		visibleSizes = new Set(visibleSizes); // Trigger reactivity
	}

	function toggleEdges() {
		showEdges = !showEdges;
		if (showEdges) {
			// Add edges for all visible sizes
			for (const size of visibleSizes) {
				addEdgesForSize(size);
			}
		} else {
			// Remove all edges
			for (const size of edgeGroups.keys()) {
				removeEdgesForSize(size);
			}
		}
	}

	// Complement relationship: n-note chord at (x,y,z) has (12-n)-note complement at (4-x, 4-y, 4-z)
	// This is because each diminished axis has 4 notes, and the complement uses the remaining notes

	function rebuildComplements() {
		if (!showComplements || !scene) return;
		removeComplements();
		addComplements();
	}

	function addComplements() {
		if (complementGroup) removeComplements();
		complementGroup = new THREE.Group();

		const seenPairs = new Set<string>();

		// For each visible size, draw lines to complement points
		for (const size of visibleSizes) {
			const coordMap = allChordsBySize.get(size)!;
			const complementSize = 12 - size;
			const complementMap = allChordsBySize.get(complementSize);

			if (!complementMap) continue;

			const color = SIZE_COLORS[size];

			for (const [, point] of coordMap) {
				const [x, y, z] = point.coord;
				const cx = 4 - x;
				const cy = 4 - y;
				const cz = 4 - z;

				// Create a unique key for this pair (order-independent)
				const pairKey = [[x, y, z].join(','), [cx, cy, cz].join(',')].sort().join('|');
				if (seenPairs.has(pairKey)) continue;
				seenPairs.add(pairKey);

				// Check if complement point exists (it should for valid complements)
				const complementKey = [cx, cy, cz].join(',');
				if (!complementMap.has(complementKey)) continue;

				// Skip self-complements (6-note chords at (2,2,2))
				if (x === cx && y === cy && z === cz) continue;

				// Create dashed line connecting point to its complement
				const points = [
					new THREE.Vector3(x, z, y), // Remap: X→X, Z→Y, Y→Z
					new THREE.Vector3(cx, cz, cy)
				];

				const lineGeometry = new THREE.BufferGeometry().setFromPoints(points);
				const lineMaterial = new THREE.LineDashedMaterial({
					color,
					transparent: true,
					opacity: 0.5,
					dashSize: 0.15,
					gapSize: 0.1
				});
				const line = new THREE.Line(lineGeometry, lineMaterial);
				line.computeLineDistances(); // Required for dashed lines
				complementGroup.add(line);

				// Add a small sphere at the midpoint to show the "inversion center"
				const midpoint = new THREE.Vector3(
					(x + cx) / 2,
					(z + cz) / 2,
					(y + cy) / 2
				);
				const midGeometry = new THREE.SphereGeometry(0.03, 8, 8);
				const midMaterial = new THREE.MeshBasicMaterial({
					color: 0xffffff,
					transparent: true,
					opacity: 0.3
				});
				const midSphere = new THREE.Mesh(midGeometry, midMaterial);
				midSphere.position.copy(midpoint);
				complementGroup.add(midSphere);
			}
		}

		// Add the center point (2,2,2) which is the inversion center for the whole space
		const centerGeometry = new THREE.SphereGeometry(0.1, 16, 16);
		const centerMaterial = new THREE.MeshBasicMaterial({
			color: 0xffffff,
			transparent: true,
			opacity: 0.6,
			wireframe: true
		});
		const centerSphere = new THREE.Mesh(centerGeometry, centerMaterial);
		centerSphere.position.set(2, 2, 2); // (2,2,2) remapped
		complementGroup.add(centerSphere);

		// Add label for center
		const centerLabel = createTextSprite('(2,2,2) center', 0xffffff);
		centerLabel.position.set(2, 2.4, 2);
		centerLabel.scale.set(1.2, 0.3, 1);
		complementGroup.add(centerLabel);

		scene.add(complementGroup);
	}

	function removeComplements() {
		if (complementGroup) {
			scene.remove(complementGroup);
			complementGroup = null;
		}
	}

	function toggleComplements() {
		showComplements = !showComplements;
		if (showComplements) {
			addComplements();
		} else {
			removeComplements();
		}
	}

	function onCanvasClick(event: MouseEvent) {
		const rect = container.getBoundingClientRect();
		mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
		mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

		raycaster.setFromCamera(mouse, camera);

		const allSpheres: THREE.Object3D[] = [];
		for (const group of sizeGroups.values()) {
			group.traverse((obj: THREE.Object3D) => {
				if (obj instanceof THREE.Mesh) {
					allSpheres.push(obj);
				}
			});
		}

		const intersects = raycaster.intersectObjects(allSpheres);
		if (intersects.length > 0) {
			const data = pointToData.get(intersects[0].object);
			if (data) {
				selectedPoint = data.point;
			}
		} else {
			selectedPoint = null;
		}
	}

	onMount(() => {
		scene = new THREE.Scene();
		scene.background = new THREE.Color(0x1a1a2e);

		camera = new THREE.PerspectiveCamera(
			75,
			container.clientWidth / container.clientHeight,
			0.1,
			1000
		);
		camera.position.set(6, 6, 6);

		raycaster = new THREE.Raycaster();
		mouse = new THREE.Vector2();

		const renderer = new THREE.WebGLRenderer({ antialias: true });
		renderer.setSize(container.clientWidth, container.clientHeight);
		renderer.setPixelRatio(window.devicePixelRatio);
		container.appendChild(renderer.domElement);

		const controls = new OrbitControls(camera, renderer.domElement);
		controls.enableDamping = true;
		controls.dampingFactor = 0.05;

		// Create axes
		// Mapping: X (Dim C) → 3D X (right), Z (Dim E) → 3D Y (up), Y (Dim D) → 3D Z (depth)
		const axisLength = 5;
		const axisRadius = 0.02;

		// X axis (red) - Dim C - goes RIGHT
		const xGeometry = new THREE.CylinderGeometry(axisRadius, axisRadius, axisLength, 16);
		const xMaterial = new THREE.MeshBasicMaterial({ color: 0xff4444 });
		const xAxis = new THREE.Mesh(xGeometry, xMaterial);
		xAxis.rotation.z = -Math.PI / 2;
		xAxis.position.x = axisLength / 2;
		scene.add(xAxis);

		// Z axis (blue) - Dim E - goes UP
		const zGeometry = new THREE.CylinderGeometry(axisRadius, axisRadius, axisLength, 16);
		const zMaterial = new THREE.MeshBasicMaterial({ color: 0x4444ff });
		const zAxis = new THREE.Mesh(zGeometry, zMaterial);
		zAxis.position.y = axisLength / 2;
		scene.add(zAxis);

		// Y axis (green) - Dim D - goes into DEPTH
		const yGeometry = new THREE.CylinderGeometry(axisRadius, axisRadius, axisLength, 16);
		const yMaterial = new THREE.MeshBasicMaterial({ color: 0x44ff44 });
		const yAxis = new THREE.Mesh(yGeometry, yMaterial);
		yAxis.rotation.x = Math.PI / 2;
		yAxis.position.z = axisLength / 2;
		scene.add(yAxis);

		// Arrow heads
		const coneRadius = 0.1;
		const coneHeight = 0.25;

		const xCone = new THREE.Mesh(new THREE.ConeGeometry(coneRadius, coneHeight, 16), xMaterial);
		xCone.rotation.z = -Math.PI / 2;
		xCone.position.x = axisLength + coneHeight / 2;
		scene.add(xCone);

		const zCone = new THREE.Mesh(new THREE.ConeGeometry(coneRadius, coneHeight, 16), zMaterial);
		zCone.position.y = axisLength + coneHeight / 2;
		scene.add(zCone);

		const yCone = new THREE.Mesh(new THREE.ConeGeometry(coneRadius, coneHeight, 16), yMaterial);
		yCone.rotation.x = Math.PI / 2;
		yCone.position.z = axisLength + coneHeight / 2;
		scene.add(yCone);

		// Axis labels
		const xLabel = createTextSprite('X: C Eb Gb A', 0xff4444);
		xLabel.position.set(axisLength + 1, 0, 0);
		xLabel.scale.set(1.8, 0.45, 1);
		scene.add(xLabel);

		const zLabel = createTextSprite('Z: E G Bb Db', 0x4444ff);
		zLabel.position.set(0, axisLength + 0.6, 0);
		zLabel.scale.set(1.8, 0.45, 1);
		scene.add(zLabel);

		const yLabel = createTextSprite('Y: D F Ab B', 0x44ff44);
		yLabel.position.set(0, 0, axisLength + 1);
		yLabel.scale.set(1.8, 0.45, 1);
		scene.add(yLabel);

		// Origin
		const originGeometry = new THREE.SphereGeometry(0.08, 16, 16);
		const originMaterial = new THREE.MeshBasicMaterial({ color: 0xffffff });
		const origin = new THREE.Mesh(originGeometry, originMaterial);
		scene.add(origin);

		// Add initial visible sizes
		for (const size of visibleSizes) {
			addSizeGroupToScene(size);
		}

		// Click handler
		container.addEventListener('click', onCanvasClick);

		function animate() {
			requestAnimationFrame(animate);
			controls.update();
			renderer.render(scene, camera);
		}
		animate();

		const handleResize = () => {
			camera.aspect = container.clientWidth / container.clientHeight;
			camera.updateProjectionMatrix();
			renderer.setSize(container.clientWidth, container.clientHeight);
		};
		window.addEventListener('resize', handleResize);

		return () => {
			window.removeEventListener('resize', handleResize);
			container.removeEventListener('click', onCanvasClick);
			renderer.dispose();
			container.removeChild(renderer.domElement);
		};
	});

	// Reactive effect: rebuild complements when visible sizes change
	$effect(() => {
		// Access visibleSizes to create dependency
		const sizes = Array.from(visibleSizes);
		// Only rebuild if complements are enabled and scene exists
		if (showComplements && scene && sizes.length >= 0) {
			rebuildComplements();
		}
	});

	// Stats for each size
	function getSizeStats(size: number) {
		const coordMap = allChordsBySize.get(size)!;
		const totalChords = combinations(ALL_NOTES, size).length;
		const uniquePoints = coordMap.size;
		return { totalChords, uniquePoints };
	}
</script>

<div class="app">
	<div bind:this={container} class="canvas-container"></div>

	<div class="panel">
		<h2>Function Cube</h2>
		<p class="subtitle">All Note Combinations</p>

		<div class="legend">
			<div class="legend-item">
				<span class="axis-dot x"></span>
				<span>X: Dim C (C, Eb, Gb, A) → right</span>
			</div>
			<div class="legend-item">
				<span class="axis-dot z"></span>
				<span>Z: Dim E (E, G, Bb, Db) → up</span>
			</div>
			<div class="legend-item">
				<span class="axis-dot y"></span>
				<span>Y: Dim D (D, F, Ab, B) → depth</span>
			</div>
		</div>

		<div class="size-selector">
			<h3>Chord Size (# of notes)</h3>
			<div class="size-buttons">
				{#each Array.from({ length: 12 }, (_, i) => i + 1) as size}
					{@const stats = getSizeStats(size)}
					<button
						class="size-btn"
						class:active={visibleSizes.has(size)}
						style="--size-color: #{SIZE_COLORS[size].toString(16).padStart(6, '0')}"
						onclick={() => toggleSize(size)}
					>
						<span class="size-num">{size}</span>
						<span class="size-stats">{stats.totalChords} → {stats.uniquePoints}pts</span>
					</button>
				{/each}
			</div>
		</div>

		<div class="edge-toggle">
			<button class="toggle-btn" class:active={showEdges} onclick={toggleEdges}>
				{showEdges ? 'Hide' : 'Show'} Semitone Edges
			</button>
			<p class="edge-info">Edges connect chords that differ by one semitone</p>
		</div>

		<div class="complement-toggle">
			<button class="toggle-btn complement" class:active={showComplements} onclick={toggleComplements}>
				{showComplements ? 'Hide' : 'Show'} Complement Links
			</button>
			<p class="edge-info">
				Connects chords to their complements.
				n-note chord at (x,y,z) ↔ (12-n)-note complement at (4-x, 4-y, 4-z).
			</p>
		</div>

		{#if selectedPoint}
			{@const groupedByType = selectedPoint.chords.reduce(
				(acc, chord) => {
					const type = chord.chordType || 'Unknown';
					if (!acc[type]) acc[type] = [];
					acc[type].push(chord);
					return acc;
				},
				{} as Record<string, ChordData[]>
			)}
			<div class="selected-info">
				<h3>Point ({selectedPoint.coord.join(', ')})</h3>
				<p class="chord-count">{selectedPoint.count} chord(s) at this location</p>
				<div class="chord-types-summary">
					{#each Object.entries(groupedByType) as [type, chords]}
						<span class="type-badge">{type}: {chords.length}</span>
					{/each}
				</div>
				<div class="chord-list">
					{#each Object.entries(groupedByType) as [type, chords]}
						<div class="type-group">
							<div class="type-header">{type}</div>
							{#each chords.slice(0, 8) as chord}
								<div class="chord-item">
									<span class="chord-notes">{chord.notes.join(' ')}</span>
									<span class="chord-props">
										<span class="prop-pattern" title="Interval pattern (axis transitions)">{chord.intervalPattern.arrows || '—'}</span>
										<span class="prop-45" class:has-interval={chord.fourthFifth !== 'none'} title="Perfect 4th/5th">
											{#if chord.fourthFifth === 'both'}45
											{:else if chord.fourthFifth === 'P4'}4
											{:else if chord.fourthFifth === 'P5'}5
											{:else}—{/if}
										</span>
										<span class="prop-27" class:has-interval={chord.semitoneMaj7 !== 'none'} title="minor 2nd / Major 7th">
											{#if chord.semitoneMaj7 === 'both'}s7
											{:else if chord.semitoneMaj7 === 'm2'}s
											{:else if chord.semitoneMaj7 === 'M7'}7
											{:else}—{/if}
										</span>
									</span>
								</div>
							{/each}
							{#if chords.length > 8}
								<div class="chord-item more">+{chords.length - 8} more {type}</div>
							{/if}
						</div>
					{/each}
				</div>
				<div class="prop-legend">
					<div class="prop-legend-item">
						<strong>Arrows:</strong> ↑same axis, ←forward, →backward
					</div>
					<div class="prop-legend-item">
						<strong>4/5:</strong> Perfect 4th/5th present
					</div>
					<div class="prop-legend-item">
						<strong>s/7:</strong> semitone (m2) / Major 7th present
					</div>
				</div>
			</div>
		{/if}

		<div class="info">
			<p>Click spheres to see chords</p>
			<p>Sphere size = # of chords</p>
			<p>Drag to rotate, scroll to zoom</p>
		</div>

		<div class="bitmask-section">
			<button class="toggle-btn" class:active={showBitmaskDemo} onclick={() => showBitmaskDemo = !showBitmaskDemo}>
				{showBitmaskDemo ? 'Hide' : 'Show'} Bitmask Demo
			</button>

			{#if showBitmaskDemo}
				{@const transposed = demoChord.transpose(demoTranspose)}
				{@const comp = demoChord.complement()}
				{@const lookup = getLookup()}
				{@const entry = lookup.get(demoChord.bits)}
				<div class="bitmask-demo">
					<h4>u16 Bitmask Encoding</h4>

					<div class="demo-chord">
						<div class="demo-label">Current Chord</div>
						<div class="demo-notes">{demoChord.notes.join(' ')}</div>
						<div class="demo-binary">
							<span class="binary-label">bits:</span>
							<span class="binary-value">{demoChord.bits.toString(2).padStart(12, '0')}</span>
						</div>
						<div class="demo-decimal">
							<span class="decimal-label">int:</span>
							<span class="decimal-value">{demoChord.bits}</span>
						</div>
						<div class="demo-cog">
							<span class="cog-label">CoG:</span>
							<span class="cog-value">({demoChord.cog.join(', ')})</span>
						</div>
						{#if entry.name}
							<div class="demo-type">
								<span class="type-label">type:</span>
								<span class="type-value">{entry.name}</span>
							</div>
						{/if}
					</div>

					<div class="note-picker">
						<div class="picker-label">Toggle Notes:</div>
						<div class="note-buttons">
							{#each ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B'] as note}
								{@const isActive = demoChord.bits & (1 << SEMI[note as NoteName])}
								<button
									class="note-btn"
									class:active={isActive}
									onclick={() => {
										const newBits = demoChord.bits ^ (1 << SEMI[note as NoteName]);
										demoChord = new Chord(newBits);
									}}
								>
									{note}
								</button>
							{/each}
						</div>
					</div>

					<div class="transpose-section">
						<div class="picker-label">Transpose: {demoTranspose > 0 ? '+' : ''}{demoTranspose}</div>
						<input
							type="range"
							min="-11"
							max="11"
							bind:value={demoTranspose}
						/>
						{#if demoTranspose !== 0}
							<div class="transposed-result">
								<span class="result-arrow">→</span>
								<span class="result-notes">{transposed.notes.join(' ')}</span>
								<span class="result-bits">({transposed.bits})</span>
							</div>
						{/if}
					</div>

					<div class="complement-section">
						<div class="picker-label">Complement (remaining notes):</div>
						<div class="complement-result">
							<span class="comp-notes">{comp.notes.join(' ')}</span>
							<span class="comp-cog">CoG: ({comp.cog.join(', ')})</span>
						</div>
						<div class="complement-math">
							({demoChord.cog.join(', ')}) + ({comp.cog.join(', ')}) = (4, 4, 4)
						</div>
					</div>

					<div class="interval-section">
						<div class="picker-label">Interval Vector:</div>
						<div class="interval-vector">
							{#each demoChord.intervalVector as count, i}
								<div class="iv-item">
									<span class="iv-class">{['m2/M7', 'M2/m7', 'm3/M6', 'M3/m6', 'P4/P5', 'TT'][i]}</span>
									<span class="iv-count">{count}</span>
								</div>
							{/each}
						</div>
					</div>

					<div class="quick-chords">
						<div class="picker-label">Quick Load:</div>
						<div class="quick-buttons">
							<button onclick={() => demoChord = Chord.fromNotes(['C', 'E', 'G'])}>C</button>
							<button onclick={() => demoChord = Chord.fromNotes(['C', 'Eb', 'G'])}>Cm</button>
							<button onclick={() => demoChord = Chord.fromNotes(['C', 'E', 'G', 'B'])}>Cmaj7</button>
							<button onclick={() => demoChord = Chord.fromNotes(['C', 'Eb', 'Gb', 'A'])}>Cdim7</button>
							<button onclick={() => demoChord = Chord.fromNotes(['C', 'E', 'Ab'])}>Caug</button>
						</div>
					</div>
				</div>
			{/if}
		</div>
	</div>
</div>

<style>
	.app {
		display: flex;
		width: 100vw;
		height: 100vh;
		overflow: hidden;
	}

	.canvas-container {
		flex: 1;
		height: 100%;
	}

	.panel {
		width: 340px;
		background: linear-gradient(180deg, #1a1a2e 0%, #16213e 100%);
		padding: 20px;
		overflow-y: auto;
		border-left: 1px solid #333;
		color: #fff;
	}

	h2 {
		margin: 0 0 5px 0;
		font-size: 1.5rem;
		color: #fff;
	}

	.subtitle {
		margin: 0 0 20px 0;
		color: #888;
		font-size: 0.9rem;
	}

	.legend {
		background: rgba(255, 255, 255, 0.05);
		border-radius: 8px;
		padding: 12px;
		margin-bottom: 20px;
	}

	.legend-item {
		display: flex;
		align-items: center;
		gap: 10px;
		margin: 8px 0;
		font-size: 0.85rem;
	}

	.axis-dot {
		width: 12px;
		height: 12px;
		border-radius: 50%;
	}

	.axis-dot.x {
		background: #ff4444;
	}
	.axis-dot.y {
		background: #44ff44;
	}
	.axis-dot.z {
		background: #4444ff;
	}

	h3 {
		font-size: 0.9rem;
		color: #aaa;
		margin: 0 0 10px 0;
		text-transform: uppercase;
		letter-spacing: 1px;
	}

	.size-selector {
		margin-bottom: 16px;
	}

	.edge-toggle {
		margin-bottom: 16px;
		padding-bottom: 16px;
		border-bottom: 1px solid #333;
	}

	.toggle-btn {
		width: 100%;
		padding: 10px;
		border: 2px solid #888;
		background: transparent;
		color: #888;
		border-radius: 6px;
		cursor: pointer;
		font-family: inherit;
		font-size: 0.9rem;
		transition: all 0.2s;
	}

	.toggle-btn:hover {
		background: rgba(255, 255, 255, 0.1);
	}

	.toggle-btn.active {
		border-color: #ffd93d;
		color: #ffd93d;
		background: rgba(255, 217, 61, 0.1);
	}

	.toggle-btn.complement {
		border-color: #e91e63;
		color: #e91e63;
	}

	.toggle-btn.complement:hover {
		background: rgba(233, 30, 99, 0.1);
	}

	.toggle-btn.complement.active {
		border-color: #e91e63;
		color: #e91e63;
		background: rgba(233, 30, 99, 0.2);
	}

	.complement-toggle {
		margin-bottom: 16px;
		padding-bottom: 16px;
		border-bottom: 1px solid #333;
	}

	.edge-info {
		font-size: 0.75rem;
		color: #666;
		margin: 8px 0 0 0;
		text-align: center;
	}

	.size-buttons {
		display: grid;
		grid-template-columns: repeat(4, 1fr);
		gap: 6px;
	}

	.size-btn {
		display: flex;
		flex-direction: column;
		align-items: center;
		padding: 8px 4px;
		border: 2px solid var(--size-color);
		background: transparent;
		color: var(--size-color);
		border-radius: 6px;
		cursor: pointer;
		transition: all 0.2s;
		font-family: inherit;
	}

	.size-btn:hover {
		background: rgba(255, 255, 255, 0.1);
	}

	.size-btn.active {
		background: var(--size-color);
		color: #1a1a2e;
	}

	.size-num {
		font-size: 1.1rem;
		font-weight: bold;
	}

	.size-stats {
		font-size: 0.6rem;
		opacity: 0.8;
	}

	.selected-info {
		background: rgba(255, 255, 255, 0.08);
		border-radius: 8px;
		padding: 12px;
		margin-bottom: 16px;
	}

	.selected-info h3 {
		color: #fff;
		margin-bottom: 8px;
	}

	.chord-count {
		color: #888;
		font-size: 0.85rem;
		margin: 0 0 10px 0;
	}

	.chord-list {
		max-height: 200px;
		overflow-y: auto;
	}

	.chord-item {
		font-family: monospace;
		font-size: 0.8rem;
		padding: 4px 8px;
		background: rgba(0, 0, 0, 0.3);
		border-radius: 4px;
		margin: 4px 0;
	}

	.chord-item.more {
		color: #888;
		font-style: italic;
	}

	.chord-item {
		display: flex;
		justify-content: space-between;
		align-items: center;
	}

	.chord-notes {
		flex: 1;
	}

	.chord-props {
		display: flex;
		gap: 6px;
		font-size: 0.7rem;
		color: #666;
	}

	.prop-pattern {
		color: #9b9;
		min-width: 32px;
		text-align: center;
	}

	.prop-45,
	.prop-27 {
		min-width: 18px;
		text-align: center;
		padding: 1px 3px;
		border-radius: 3px;
		background: rgba(0, 0, 0, 0.3);
	}

	.prop-45.has-interval {
		color: #4d96ff;
		background: rgba(77, 150, 255, 0.2);
	}

	.prop-27.has-interval {
		color: #ff6b6b;
		background: rgba(255, 107, 107, 0.2);
	}

	.prop-legend {
		margin-top: 12px;
		padding-top: 10px;
		border-top: 1px solid #333;
		font-size: 0.7rem;
		color: #666;
	}

	.prop-legend-item {
		margin: 4px 0;
	}

	.prop-legend-item strong {
		color: #888;
	}

	.chord-types-summary {
		display: flex;
		flex-wrap: wrap;
		gap: 6px;
		margin-bottom: 12px;
	}

	.type-badge {
		font-size: 0.7rem;
		padding: 2px 8px;
		background: rgba(255, 255, 255, 0.15);
		border-radius: 10px;
		color: #ccc;
	}

	.type-group {
		margin-bottom: 10px;
	}

	.type-header {
		font-size: 0.75rem;
		font-weight: bold;
		color: #ffd93d;
		margin-bottom: 4px;
		text-transform: uppercase;
	}

	.info {
		margin-top: 20px;
		padding-top: 16px;
		border-top: 1px solid #333;
		font-size: 0.8rem;
		color: #666;
	}

	.info p {
		margin: 4px 0;
	}

	:global(body) {
		margin: 0;
		padding: 0;
		font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
	}

	/* Bitmask Demo Styles */
	.bitmask-section {
		margin-top: 16px;
		padding-top: 16px;
		border-top: 1px solid #333;
	}

	.bitmask-demo {
		margin-top: 12px;
		background: rgba(0, 0, 0, 0.3);
		border-radius: 8px;
		padding: 12px;
	}

	.bitmask-demo h4 {
		margin: 0 0 12px 0;
		color: #4d96ff;
		font-size: 0.9rem;
	}

	.demo-chord {
		background: rgba(255, 255, 255, 0.05);
		border-radius: 6px;
		padding: 10px;
		margin-bottom: 12px;
	}

	.demo-label {
		font-size: 0.7rem;
		color: #888;
		text-transform: uppercase;
		margin-bottom: 4px;
	}

	.demo-notes {
		font-size: 1.2rem;
		font-weight: bold;
		color: #fff;
		font-family: monospace;
		margin-bottom: 8px;
	}

	.demo-binary,
	.demo-decimal,
	.demo-cog,
	.demo-type {
		display: flex;
		gap: 8px;
		font-size: 0.8rem;
		margin: 4px 0;
	}

	.binary-label,
	.decimal-label,
	.cog-label,
	.type-label {
		color: #666;
		min-width: 40px;
	}

	.binary-value {
		font-family: monospace;
		color: #6bcb77;
		letter-spacing: 1px;
	}

	.decimal-value {
		font-family: monospace;
		color: #ffd93d;
	}

	.cog-value {
		font-family: monospace;
		color: #4d96ff;
	}

	.type-value {
		color: #e91e63;
		font-weight: bold;
	}

	.note-picker,
	.transpose-section,
	.complement-section,
	.interval-section,
	.quick-chords {
		margin-top: 12px;
	}

	.picker-label {
		font-size: 0.75rem;
		color: #888;
		margin-bottom: 6px;
	}

	.note-buttons {
		display: grid;
		grid-template-columns: repeat(6, 1fr);
		gap: 4px;
	}

	.note-btn {
		padding: 6px 4px;
		border: 1px solid #444;
		background: transparent;
		color: #888;
		border-radius: 4px;
		cursor: pointer;
		font-family: monospace;
		font-size: 0.75rem;
		transition: all 0.15s;
	}

	.note-btn:hover {
		background: rgba(255, 255, 255, 0.1);
	}

	.note-btn.active {
		background: #4d96ff;
		border-color: #4d96ff;
		color: #fff;
	}

	.transpose-section input[type="range"] {
		width: 100%;
		margin: 4px 0;
	}

	.transposed-result {
		display: flex;
		align-items: center;
		gap: 8px;
		font-size: 0.85rem;
		margin-top: 6px;
	}

	.result-arrow {
		color: #666;
	}

	.result-notes {
		font-family: monospace;
		color: #6bcb77;
	}

	.result-bits {
		font-family: monospace;
		color: #666;
		font-size: 0.75rem;
	}

	.complement-result {
		display: flex;
		justify-content: space-between;
		align-items: center;
		background: rgba(233, 30, 99, 0.1);
		border: 1px solid rgba(233, 30, 99, 0.3);
		border-radius: 4px;
		padding: 6px 8px;
	}

	.comp-notes {
		font-family: monospace;
		color: #e91e63;
		font-size: 0.85rem;
	}

	.comp-cog {
		font-family: monospace;
		color: #888;
		font-size: 0.75rem;
	}

	.complement-math {
		font-family: monospace;
		font-size: 0.7rem;
		color: #666;
		text-align: center;
		margin-top: 6px;
	}

	.interval-vector {
		display: grid;
		grid-template-columns: repeat(3, 1fr);
		gap: 4px;
	}

	.iv-item {
		display: flex;
		justify-content: space-between;
		background: rgba(255, 255, 255, 0.05);
		border-radius: 4px;
		padding: 4px 6px;
		font-size: 0.7rem;
	}

	.iv-class {
		color: #888;
	}

	.iv-count {
		font-family: monospace;
		color: #ffd93d;
		font-weight: bold;
	}

	.quick-buttons {
		display: flex;
		flex-wrap: wrap;
		gap: 4px;
	}

	.quick-buttons button {
		padding: 4px 10px;
		border: 1px solid #444;
		background: transparent;
		color: #888;
		border-radius: 4px;
		cursor: pointer;
		font-size: 0.75rem;
		transition: all 0.15s;
	}

	.quick-buttons button:hover {
		background: rgba(255, 255, 255, 0.1);
		color: #fff;
	}
</style>
