// ═══════════════════════════════════════════════════════════════════
//  ΔΥΑΔΙΚΟΤΗΤΑ
//  Ηλεκτροακουστική σύνθεση σε τρία μέρη για το κύμα-σωματίδιο
//
//  Α · Η ΜΕΤΑΜΟΡΦΩΣΗ — ένα ηχητικό αντικείμενο ταξιδεύει από σωματίδιο
//      σε κύμα, με στάση στον «καθρέφτη» όπου είναι και τα δύο
//  Β · ΤΟ ΕΡΓΑΣΤΗΡΙΟ — το πείραμα της διπλής σχισμής: μέτρηση, συμβολή,
//      κατάρρευση, ελευθερία· το κυματικό νήμα αντιδρά στην παρατήρηση
//  Γ · ΤΟ ΦΑΣΜΑ — ηχοποίηση της σειράς Balmer του υδρογόνου
//      (πηγή δεδομένων: NIST Atomic Spectra Database — Kramida,
//      Ralchenko, Reader & NIST ASD Team) σε Λα μινόρε πεντατονική,
//      με τους συντελεστές Einstein A ως πιθανότητες· τέλος: ιονισμός
//
//  ΑΠΑΙΤΗΣΕΙΣ: SuperCollider 3.x — τίποτα άλλο.
//  ΕΚΤΕΛΕΣΗ: κέρσορας οπουδήποτε ΜΕΣΑ στο μπλοκ, Ctrl+Enter
//            (Windows/Linux) ή Cmd+Enter (macOS). Τίποτα άλλο.
//  ΔΙΑΚΟΠΗ: Ctrl+.  (Cmd+.)
//
//  Δημιουργήθηκε σε διάλογο με το Claude (Anthropic) — η πλήρης
//  τεκμηρίωση της διαδικασίας βρίσκεται στον φάκελο queries/
// ═══════════════════════════════════════════════════════════════════
(
s.waitForBoot {

	// ═══════════ ΧΡΟΝΟΣ — όλες οι αποφάσεις διάρκειας εδώ ═══════════
	~z = 0.6;   // ΚΛΙΜΑΚΑ ΧΡΟΝΟΥ ΟΛΟΥ ΤΟΥ ΕΡΓΟΥ
	            //   0.6 → ~10 λεπτά (τρέχουσα εκδοχή)
	            //   1.0 → ~17 λεπτά (long-form: Α=4', φάσεις Β=2')
	            //   0.2 → ~3.5 λεπτά (γρήγορη πρόβα δομής)
	// Αν το ζητούμενο ήταν ΜΟΝΟ το Γ ≈ 10 λεπτά:
	//   βάλε ~z = 1 και ~durC = 600 * ~z

	~durA  = 240 * ~z;   // Μέρος Α · η μεταμόρφωση
	~gefAB = 6   * ~z;   // γέφυρα Α→Β: σιωπή — η πρώτη μέτρηση τη σπάει
	~durP1 = 120 * ~z;   // Β/1 · μέτρηση
	~durTr = 60  * ~z;   // Β/μετάβαση Gabor
	~durP2 = 120 * ~z;   // Β/2 · δύο σχισμές
	~durP3 = 120 * ~z;   // Β/3 · κατάρρευση
	~durP4 = 120 * ~z;   // Β/4 · ελευθερία
	~gefBC = 7   * ~z;   // ανάσα Β→Γ
	~durC  = 150 * ~z;   // Μέρος Γ (χωρίς την ουρά του ιονισμού)

	// ═══════════ ΥΛΙΚΟ ═══════════
	~fund    = 110;      // θεμέλιος του Β — οκτάβα κάτω από το Λα (220) των Α/Γ
	~detLo   = 660;      // ο ανιχνευτής του Β κληρώνει συχνότητα
	~detHi   = 760;      //   μέσα σε αυτό το παράθυρο
	~temper  = 0.4;      // Γ: τιθάσευση των πιθανοτήτων Einstein (1 = ωμή φύση)
	~isopedo = 0.05;     // Γ: έδαφος 220/220.5 — ο πρώτος ήχος του project, επιστρέφει

	// Δεδομένα: NIST ASD — υδρογόνο, σειρά Balmer (μεταβάσεις n→2)
	~nm   = [656.28, 486.13, 434.05, 410.17, 397.01, 388.90]; // μήκη κύματος (nm)
	~aki  = [44.1, 8.42, 2.53, 0.97, 0.44, 0.22];             // Einstein A (×10⁶ s⁻¹)
	~nHi  = [3, 4, 5, 6, 7, 8];                                // αρχική στάθμη κάθε άλματος
	~wq   = (~aki ** ~temper).normalizeSum;                    // τα «ζάρια» του Γ/2
	~penta = [0, 3, 5, 7, 10];                                 // Λα μινόρε πεντατονική
	~n2f  = { |n| 220 * (2 ** (~penta.wrapAt(n - 2) / 12)) * (2 ** ((n - 2).div(5))) };
	~orio = 299792458 / 364.6e-9 / (2 ** 40); // όριο Balmer ÷ 2^40 = 747.8 Hz — ο ιονισμός

	// ═══════════ ΤΑ ΟΡΓΑΝΑ ═══════════

	// Μέρος Α ολόκληρο ως ένας synth: αραιές ανιχνεύσεις (Poisson ~5"),
	// στάση στον καθρέφτη, τρεις αναδύσεις με δύο ανεπιτήρητες σιωπές
	SynthDef(\qm_partA, { |out = 0, dur = 240|
		var araiotita = 0.2, siopi = 0;
		var morph, amp, densP, durP, trigP, particle;
		var densW, durW, trigW, wave, beat, mixWave, sig;
		morph = EnvGen.kr(Env([0, 0.5, 0.5, 1], [0.4, 0.2, 0.4] * dur, \sin));
		amp = EnvGen.kr(Env(
			[0, 1, 1, siopi, siopi, 1, 1, siopi, siopi, 1, 1, 0],
			[0.06, 0.20, 0.04, 0.03, 0.04, 0.15, 0.04, 0.03, 0.04, 0.27, 0.10] * dur,
			\sin
		), doneAction: 2);
		densP = morph.linexp(0, 1, araiotita, 45);
		durP  = morph.linlin(0, 1, 0.09, 0.03);
		trigP = Dust.ar(densP);
		particle = GrainSin.ar(2, trigP, durP, 1880, LFNoise0.ar(densP));
		densW = morph.linexp(0, 1, 8, 90);
		durW  = morph.linlin(0, 1, 0.03, 0.12);
		trigW = Impulse.ar(densW);
		wave  = GrainSin.ar(2, trigW, durW, 220, 0);
		beat  = Mix(SinOsc.ar([20, 20.5], 0, 0.5)) ! 2;
		mixWave = XFade2.ar(wave, beat, morph.linlin(0.65, 1, -1, 1).clip(-1, 1));
		sig = XFade2.ar(particle, mixWave, morph.linlin(0, 0.85, -1, 1).clip(-1, 1));
		Out.ar(out, sig * 0.2 * amp);
	}).add;

	// Το κυματικό νήμα του Β — ζαρώνει όταν το μετρούν, ανθίζει όταν όχι
	SynthDef(\qm_waveState, {
		|out = 0, fund = 110, amp = 0.22, gate = 1, glide = 6, rel = 8|
		var f = Lag.kr(fund, glide);
		var a = Lag.kr(amp, glide);
		var ratios  = [1, 1.5, 2, 2.5, 3, 4, 5, 6];
		var weights = [1, 0.7, 0.55, 0.4, 0.3, 0.2, 0.14, 0.09];
		var partials = 8.collect({ |i|
			var freq   = f * ratios[i];
			var detune = LFNoise1.kr(0.05 + (i * 0.02)).range(-0.8, 0.8);
			var pan    = LFNoise1.kr(0.07 + (i * 0.015)).range(-1, 1);
			Pan2.ar(SinOsc.ar(freq + detune) * weights[i], pan)
		}).sum;
		var shimmer = LFNoise2.kr(0.15).range(0.75, 1.0);
		var env     = EnvGen.kr(Env.asr(4.0, 1, rel), gate, doneAction: Done.freeSelf);
		Out.ar(out, partials * a * shimmer * env);
	}).add;

	SynthDef(\qm_waveCloud, {
		|out = 0, fund = 55, amp = 0.1, gate = 1|
		var lfo = LFNoise2.kr(0.08).range(0.98, 1.02);
		var sig = Mix.fill(5, { |i|
			var f = fund * (i + 1) * lfo * LFNoise1.kr(0.03).range(0.995, 1.005);
			var p = LFNoise1.kr(0.04 + (i * 0.01)).range(-1, 1);
			Pan2.ar(SinOsc.ar(f) * 0.2, p)
		});
		var env = EnvGen.kr(Env.asr(6.0, 1, 10.0), gate, doneAction: Done.freeSelf);
		Out.ar(out, sig * amp * env);
	}).add;

	// Το «κλικ» του ανιχνευτή
	SynthDef(\qm_particle, {
		|out = 0, freq = 440, pan = 0, amp = 0.7, decay = 0.5|
		var env = EnvGen.kr(
			Env.new([0, 1, 0.2, 0], [0.003, decay * 0.25, decay * 0.75]),
			doneAction: Done.freeSelf
		);
		Out.ar(out, Pan2.ar((Impulse.ar(0) * 0.3 + (SinOsc.ar(freq) * env * 0.85)) * amp, pan));
	}).add;

	SynthDef(\qm_reform, {
		|out = 0, fund = 110, amp = 0.12, gate = 1|
		var sig = Mix.fill(6, { |i|
			var freq = fund * (i + 1) * LFNoise1.kr(0.3).range(0.99, 1.01);
			Pan2.ar(SinOsc.ar(freq), LFNoise1.kr(0.5).range(-1, 1)) * (1 / (i + 1))
		});
		Out.ar(out, sig * EnvGen.kr(Env.asr(0.5, 1, 2.5), gate, doneAction: Done.freeSelf) * amp);
	}).add;

	// Η αρχή της αβεβαιότητας (όριο Gabor): βέβαιη συχνότητα = μακρύς ήχος
	SynthDef(\qm_uncertainty, {
		|out = 0, certainty = 0.5, freq = 330, amp = 0.4|
		var dur = certainty.linexp(0, 1, 0.02, 8.0);
		var rq  = certainty.linexp(0, 1, 1.8, 0.003);
		var env = EnvGen.kr(Env.sine(dur), doneAction: Done.freeSelf);
		Out.ar(out, Pan2.ar(BPF.ar(WhiteNoise.ar, freq, rq) * env * amp, 0));
	}).add;

	// Οι κροσσοί: κάθε αρμονική έχει τους δικούς της, όπως κάθε μήκος κύματος
	SynthDef(\qm_interferenceRich, {
		|out = 0, freqBase = 110, phaseOffset = 0, amp = 0.22|
		var layers = [1, 2, 3].collect({ |mult|
			var f     = freqBase * mult;
			var slit1 = SinOsc.ar(f, 0);
			var slit2 = SinOsc.ar(f, phaseOffset * mult);
			var w     = 1 / mult;
			[(slit1 + slit2) * w * 0.5, (slit1 - slit2) * w * 0.5]
		});
		var sums  = layers.collect({ |l| l[0] }).sum;
		var diffs = layers.collect({ |l| l[1] }).sum;
		Out.ar(out, [(sums + (diffs * 0.5)) * amp, (sums - (diffs * 0.5)) * amp]);
	}).add;

	// Γ: το «φωτόνιο» — κάθε νότα ένα σωματιδιακό γεγονός
	SynthDef(\c_pluck, { |out = 0, freq = 440, amp = 0.3, pan = 0, rel = 1.5|
		var sig = SinOsc.ar(freq) + SinOsc.ar(freq * 2, 0, 0.25);
		Out.ar(out, Pan2.ar(sig * EnvGen.kr(Env.perc(0.004, rel), doneAction: 2) * amp, pan));
	}).add;

	// Γ: το έδαφος — η συμβολή 220/220.5, ο πρώτος ήχος του project
	SynthDef(\c_edafos, { |out = 0, amp = 0.05, gate = 1, atk = 15, rel = 12|
		var sig = SinOsc.ar([220, 220.5]);
		var env = EnvGen.kr(Env.asr(atk, 1, rel), gate, doneAction: 2);
		Out.ar(out, sig * amp * env);
	}).add;

	s.sync;

	// ═══════════ ΒΟΗΘΟΙ ═══════════
	// Poisson αναμονή με μέσο όρο mean και ταβάνι ασφαλείας cap
	~expWait  = { |mean, cap = 12| (1.0.rand.log.neg * mean).min(cap) };
	~fadeFree = { |syn, from = 0.22, time = 1.5|
		Routine({
			30.do { |i| syn.set(\amp, from * (1 - ((i + 1) / 30))); (time / 30).wait };
			syn.free;
		}).play;
	};

	~est = ~durA + ~gefAB + (4 * ~z) + ~durP1 + ~durTr + ~durP2 + ~durP3
		+ ~durP4 + (54 * ~z) + ~gefBC + (~durC * 0.85) + (25 * ~z) + 8;
	("ΔΥΑΔΙΚΟΤΗΤΑ — όργανα έτοιμα. Εκτιμώμενη διάρκεια: ~"
		++ (~est / 60).round(0.1) ++ " λεπτά").postln;

	// ═══════════ Η ΠΑΡΤΙΤΟΥΡΑ ═══════════
	Routine({
		var ws, interf, sweep, cloud, reform1, reform2, edafos;
		var t, w, i, mean;

		// ───────────── ΜΕΡΟΣ Α · Η ΜΕΤΑΜΟΡΦΩΣΗ ─────────────
		"ΜΕΡΟΣ Α · η μεταμόρφωση".postln;
		Synth(\qm_partA, [\dur, ~durA]);

		// ─── γέφυρα Α→Β: σιωπή — και η πρώτη μέτρηση τη σπάει
		(~durA + ~gefAB).wait;

		// ───────────── ΜΕΡΟΣ Β · ΤΟ ΕΡΓΑΣΤΗΡΙΟ ─────────────
		"ΜΕΡΟΣ Β · το εργαστήριο".postln;
		ws = Synth(\qm_waveState, [\fund, ~fund, \amp, 0.05, \glide, 6, \rel, 14 * ~z]);
		(4 * ~z).wait;

		// Β/1 · ΜΕΤΡΗΣΗ — Geiger που σφίγγει ανεπαίσθητα (μέση παύση 4.5" → 1.8")
		"Β/1: παρατήρηση".postln;
		t = 0;
		while { t < ~durP1 } {
			Synth(\qm_particle, [\freq, rrand(~detLo, ~detHi), \pan, [-1, 1].choose,
				\decay, rrand(0.3, 0.6), \amp, 0.4]);
			w = ~expWait.(t.linlin(0, ~durP1, 4.5, 1.8) * ~z);
			w.wait;  t = t + w;
		};

		// Β/μετάβαση · Gabor — 9 σκαλιά από τον κρότο στον τόνο
		"Β/μετάβαση: Gabor".postln;
		ws.set(\amp, 0.11);
		9.do { |k|
			var c = k.linlin(0, 8, 0.03, 0.97);
			Synth(\qm_uncertainty, [\certainty, c, \freq, ~fund * 3, \amp, 0.4 + (c * 0.5)]);
			(~durTr / 9 * rrand(0.7, 1.3)).wait;
		};

		// Β/2 · ΔΥΟ ΣΧΙΣΜΕΣ — το νήμα ανθίζει, το πέτασμα περιπλανιέται
		"Β/2: δύο σχισμές".postln;
		ws.set(\glide, 40 * ~z);
		ws.set(\amp, 0.2);
		reform1 = Synth(\qm_reform, [\fund, ~fund, \amp, 0.14]);
		(2 * ~z).wait;
		interf = Synth(\qm_interferenceRich, [\freqBase, ~fund, \phaseOffset, 0, \amp, 0.22]);
		sweep = Routine({
			var phi = 0, inc = 0.007;
			inf.do {
				inc = (inc + 0.0015.rand2).clip(0.003, 0.012);
				phi = phi + inc;
				interf.set(\phaseOffset, phi);
				0.05.wait
			}
		}).play;
		(~durP2 - (2 * ~z)).wait;

		// Β/3 · ΚΑΤΑΡΡΕΥΣΗ — πυκνή μέτρηση που χάνει σιγά σιγά το ενδιαφέρον της
		"Β/3: κατάρρευση".postln;
		sweep.stop;
		~fadeFree.(interf, 0.22, 0.08);
		reform1.set(\gate, 0);
		ws.set(\glide, 0.4);
		ws.set(\amp, 0.045);
		Synth(\qm_particle, [\freq, rrand(~detLo, ~detHi), \pan, [-1, 1].choose,
			\decay, 0.7, \amp, 0.55]);
		(1 * ~z).wait;               // η ζάρωση προλαβαίνει να συμβεί...
		ws.set(\glide, ~durP3);      // ...και μετά το νήμα έρπει αόρατα προς τα πάνω
		ws.set(\amp, 0.09);
		t = 0;
		while { t < ~durP3 } {
			Synth(\qm_particle, [\freq, rrand(~detLo, ~detHi), \pan, [-1, 1].choose,
				\decay, rrand(0.2, 0.4), \amp, 0.35]);
			w = ~expWait.(t.linlin(0, ~durP3, 0.45, 3.5) * ~z);
			w.wait;  t = t + w;
		};

		// Β/4 · ΕΛΕΥΘΕΡΙΑ — μια οκτάβα κάτω
		"Β/4: ελευθερία".postln;
		ws.set(\glide, 10);
		ws.set(\amp, 0.22, \fund, ~fund * 0.5);
		reform2 = Synth(\qm_reform, [\fund, ~fund * 0.5, \amp, 0.12]);
		(1.5 * ~z).wait;
		interf = Synth(\qm_interferenceRich, [\freqBase, ~fund * 0.5, \phaseOffset, 0, \amp, 0.2]);
		sweep = Routine({
			var phi = 0, inc = 0.004;
			inf.do {
				inc = (inc + 0.001.rand2).clip(0.002, 0.007);
				phi = phi + inc;
				interf.set(\phaseOffset, phi);
				0.05.wait
			}
		}).play;
		cloud = Synth(\qm_waveCloud, [\fund, 55, \amp, 0.09]);
		(~durP4 - (1.5 * ~z)).wait;

		// Β/τέλος · ΑΝΕΠΙΤΗΡΗΤΟ — το εργαστήριο αποσύρεται, το κύμα βυθίζεται
		"Β/τέλος: ανεπιτήρητο".postln;
		sweep.stop;
		reform2.set(\gate, 0);
		cloud.set(\gate, 0);
		~fadeFree.(interf, 0.2, 8);
		(12 * ~z).wait;              // ολομόναχο στη σκηνή
		ws.set(\glide, 6);
		ws.set(\amp, 0.26);          // η τελευταία ανάσα — κανείς δεν κοιτάζει
		(8 * ~z).wait;
		ws.set(\glide, 18 * ~z);
		ws.set(\fund, 40, \amp, 0.0);
		(18 * ~z).wait;
		ws.set(\gate, 0);
		(16 * ~z).wait;

		// ─── ανάσα Β→Γ
		~gefBC.wait;

		// ───────────── ΜΕΡΟΣ Γ · ΤΟ ΦΑΣΜΑ ─────────────
		"ΜΕΡΟΣ Γ · το φάσμα".postln;
		(~isopedo > 0).if {
			edafos = Synth(\c_edafos, [\amp, ~isopedo, \atk, 15 * ~z, \rel, 12 * ~z]);
		};

		// Γ/1 · ΘΕΜΕΛΙΩΔΗΣ ΚΑΤΑΣΤΑΣΗ — μόνο το n=2, αραιό, μέσα από τη σιωπή
		"Γ/1: θεμελιώδης κατάσταση".postln;
		t = 0;
		while { t < (~durC * 0.22) } {
			Synth(\c_pluck, [\freq, 220, \amp, 0.22, \rel, 2.2]);
			w = ~expWait.(3.5 * ~z, 9);
			w.wait;  t = t + w;
		};

		// Γ/2 · ΟΙ ΜΕΤΑΒΑΣΕΙΣ — άλματα n→2 με τα τιθασευμένα ζάρια της φύσης
		"Γ/2: μεταβάσεις".postln;
		while { t < (~durC * 0.85) } {
			i = (0..5).wchoose(~wq);
			Synth(\c_pluck, [\freq, ~n2f.(~nHi[i]), \amp, 0.3, \pan, 1.0.rand2]);
			(0.14 * ~z).wait;
			Synth(\c_pluck, [\freq, 220, \amp, 0.2, \rel, 2]);   // κάθε άλμα γυρνά σπίτι
			mean = t.linlin(~durC * 0.22, ~durC * 0.85, 2.8, 1.1);
			w = ~expWait.(mean * ~z, 7);
			w.wait;  t = t + w;
		};

		// Γ/3 · ΑΝΑΡΡΙΧΗΣΗ — ΙΟΝΙΣΜΟΣ: τα ζάρια αντιστρέφονται, κερδίζουν οι σπάνιες στάθμες
		"Γ/3: αναρρίχηση — ιονισμός".postln;
		8.do { |k|
			var bias = k.linlin(0, 7, ~temper, -0.6);
			i = (0..5).wchoose((~aki ** bias).normalizeSum);
			Synth(\c_pluck, [\freq, ~n2f.(~nHi[i]), \amp, 0.3, \pan, 1.0.rand2]);
			(0.14 * ~z).wait;
			Synth(\c_pluck, [\freq, 220, \amp, 0.18, \rel, 2]);
			(rrand(1.2, 2.2) * ~z).wait;
		};
		edafos.notNil.if { edafos.set(\gate, 0) };
		(3 * ~z).wait;

		// Η ΔΡΑΠΕΤΕΥΣΗ: 747.8 Hz — η μόνη νότα εκτός κλίμακας.
		// Το ηλεκτρόνιο φεύγει από το άτομο ΚΑΙ από την ανθρώπινη κλίμακα.
		// Δεν προσγειώνεται ποτέ στο 220.
		Synth(\c_pluck, [\freq, ~orio, \amp, 0.34, \rel, 8]);
		(8 + (6 * ~z)).wait;
		"— ΤΕΛΟΣ ΕΡΓΟΥ —".postln;
	}).play;
};
)
