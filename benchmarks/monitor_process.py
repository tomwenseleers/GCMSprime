"""Optional Windows process-tree monitor; standard Python library only.

Usage: python benchmarks/monitor_process.py OUTPUT_PREFIX -- COMMAND [ARGS...]
Run R benchmark stages from the GCMSprime repository root.
"""
import ctypes as C
from ctypes import wintypes as W
from pathlib import Path
import csv,json,os,subprocess,sys,time

if os.name != "nt":raise SystemExit("This optional process monitor requires Windows.")
k=C.WinDLL('kernel32',use_last_error=True);ps=C.WinDLL('psapi',use_last_error=True)
class ProcessEntry(C.Structure):
    _fields_=[('dwSize',W.DWORD),('cntUsage',W.DWORD),('pid',W.DWORD),
        ('heap',C.c_size_t),('module',W.DWORD),('threads',W.DWORD),
        ('parent',W.DWORD),('priority',W.LONG),('flags',W.DWORD),('exe',W.WCHAR*260)]
class Memory(C.Structure):
    _fields_=[('cb',W.DWORD),('faults',W.DWORD)]+[(n,C.c_size_t) for n in
        ('peak_working','working','peak_paged','paged','peak_nonpaged','nonpaged','pagefile','peak_pagefile','private')]
class Status(C.Structure):
    _fields_=[('length',W.DWORD),('load',W.DWORD)]+[(n,C.c_ulonglong) for n in
        ('total','available','total_page','available_page','total_virtual','available_virtual','extended')]
k.CreateToolhelp32Snapshot.argtypes=[W.DWORD,W.DWORD];k.CreateToolhelp32Snapshot.restype=W.HANDLE
k.Process32FirstW.argtypes=[W.HANDLE,C.POINTER(ProcessEntry)];k.Process32FirstW.restype=W.BOOL
k.Process32NextW.argtypes=k.Process32FirstW.argtypes;k.Process32NextW.restype=W.BOOL
k.OpenProcess.argtypes=[W.DWORD,W.BOOL,W.DWORD];k.OpenProcess.restype=W.HANDLE
k.CloseHandle.argtypes=[W.HANDLE];k.CloseHandle.restype=W.BOOL
ps.GetProcessMemoryInfo.argtypes=[W.HANDLE,C.POINTER(Memory),W.DWORD];ps.GetProcessMemoryInfo.restype=W.BOOL
k.GetProcessTimes.argtypes=[W.HANDLE]+[C.POINTER(W.FILETIME)]*4;k.GetProcessTimes.restype=W.BOOL
k.GlobalMemoryStatusEx.argtypes=[C.POINTER(Status)];k.GlobalMemoryStatusEx.restype=W.BOOL
def processes():
    snapshot=k.CreateToolhelp32Snapshot(2,0)
    entry=ProcessEntry();entry.dwSize=C.sizeof(entry);out={}
    good=k.Process32FirstW(snapshot,C.byref(entry))
    while good:
        out[entry.pid]=(entry.parent,entry.threads,str(entry.exe))
        good=k.Process32NextW(snapshot,C.byref(entry))
    k.CloseHandle(snapshot)
    return out
def usage(pid):
    handle=k.OpenProcess(0x410,False,pid)
    if not handle:return None
    m=Memory();m.cb=C.sizeof(m)
    ok=ps.GetProcessMemoryInfo(handle,C.byref(m),C.sizeof(m))
    times=[W.FILETIME() for _ in range(4)];cpu=0
    if k.GetProcessTimes(handle,*[C.byref(t) for t in times]):
        cpu=sum((t.dwHighDateTime<<32)+t.dwLowDateTime for t in times[2:])/1e7
    k.CloseHandle(handle)
    return (m.working,m.private,cpu) if ok else None

import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("output",help="Output prefix for .log, -memory.csv and -memory.json")
parser.add_argument("command",nargs=argparse.REMAINDER)
args=parser.parse_args();command=args.command
if command and command[0]=="--":command=command[1:]
if not command:parser.error("Supply a command after --")
prefix=Path(args.output).resolve();prefix.parent.mkdir(parents=True,exist_ok=True)
env=os.environ.copy()
env.update(OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",OMP_NUM_THREADS="1",OMP_WAIT_POLICY="PASSIVE")
begun=time.monotonic();rows=[];cpu_max={}
with Path(str(prefix)+".log").open("w") as log:
    child=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=0x08000000)
    while True:
        active=processes();ids={child.pid};previous=set()
        while previous!=ids:
            previous=set(ids);ids.update(pid for pid,info in active.items() if info[0] in ids)
        working=private=threads=count=0
        for pid in ids:
            if pid not in active:continue
            value=usage(pid)
            if value:
                working+=value[0];private+=value[1];cpu_max[pid]=max(cpu_max.get(pid,0),value[2])
                threads+=active[pid][1];count+=1
        memory=Status();memory.length=C.sizeof(memory);k.GlobalMemoryStatusEx(C.byref(memory))
        rows.append(dict(seconds=time.monotonic()-begun,processes=count,threads=threads,
            working_GiB=working/1024**3,private_GiB=private/1024**3,
            available_GiB=memory.available/1024**3,cpu_seconds=sum(cpu_max.values())))
        if child.poll() is not None:break
        time.sleep(1)
elapsed=time.monotonic()-begun
with Path(str(prefix)+"-memory.csv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
summary=dict(exit_code=child.returncode,elapsed=elapsed,
    peak_working_GiB=max(x["working_GiB"] for x in rows),
    peak_private_GiB=max(x["private_GiB"] for x in rows),
    minimum_available_GiB=min(x["available_GiB"] for x in rows),
    maximum_processes=max(x["processes"] for x in rows),
    cpu_seconds=sum(cpu_max.values()),average_cpu_cores=sum(cpu_max.values())/elapsed)
Path(str(prefix)+"-memory.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2));raise SystemExit(child.returncode)
