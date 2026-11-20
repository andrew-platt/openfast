"""Run one shard of the CI regression tests serially, timing each with /usr/bin/time.

usage: python profile_serial.py <tests.tsv> <prior_costs.csv> <shard> <nshards> <out.csv>
Run from the build directory. Shards are balanced by LPT on the prior (historical) cost.
Writes CSV: ci_job,test,status,wall_s,maxrss_mb,also_in
"""
import csv, os, re, subprocess, sys, time

tsv, prior_f, shard, nshards, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]

jobs_of = {}
for line in open(tsv):
    job, test = line.rstrip("\n").split("\t")
    jobs_of.setdefault(test, []).append(job)
prior = {r["test"]: float(r["median_sec"]) for r in csv.DictReader(open(prior_f))}

# LPT assignment, deterministic so every shard computes the same partition
loads = [0.0] * nshards
assign = {}
for t in sorted(jobs_of, key=lambda t: (-prior.get(t, 10.0), t)):
    k = loads.index(min(loads))
    assign[t] = k
    loads[k] += prior.get(t, 10.0)
mine = [t for t in sorted(jobs_of, key=lambda t: -prior.get(t, 10.0)) if assign[t] == shard]
print(f"shard {shard}/{nshards}: {len(mine)} tests, prior load {loads[shard]:.0f}s", flush=True)

w = csv.writer(open(out, "w", newline=""))
w.writerow(["ci_job", "test", "status", "wall_s", "maxrss_mb", "also_in"])
for t in mine:
    job = jobs_of[t][0]
    env = dict(os.environ)
    if job == "rtest-FF":            # match CI: FF job uses 2 OpenMP threads
        env["OMP_NUM_THREADS"] = "2"
    tf = "/tmp/time.txt"
    t0 = time.time()
    p = subprocess.run(["/usr/bin/time", "-f", "%e %M", "-o", tf,
                        "ctest", "-R", f"^{re.escape(t)}$", "--timeout", "3600",
                        "--output-on-failure"],
                       env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    wall = time.time() - t0
    try:
        e, m = open(tf).read().split()[-2:]
        rss = int(m) / 1024
    except Exception:
        rss = -1
    status = "pass" if p.returncode == 0 else f"fail({p.returncode})"
    if p.returncode != 0:
        print(p.stdout[-3000:])
    print(f"{job:24s} {t:45s} {status:8s} {wall:8.1f}s {rss:8.0f}MB", flush=True)
    w.writerow([job, t, status, round(wall, 2), round(rss), ";".join(jobs_of[t][1:])])
