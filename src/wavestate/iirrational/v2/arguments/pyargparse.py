# -*- coding: utf-8 -*-
"""
"""

import argparse
from ...utilities.strings import padding_remove


def kwdict_argparse(ap, kwdict, groups_kw = dict()):
    groups = dict()
    groups_prio = dict()
    for gname, gdict in groups_kw.items():
        APpriority = gdict.get('APpriority', None)
        if APpriority is not None:
            groups_prio[gname] = APpriority

    for hname, hdict in kwdict.items():
        if hdict.get('APignore', False):
            continue
        group = hdict.get('APgroup', None)
        s = groups.setdefault(group, set())
        s.add(hname)

        gprio = groups_prio.get(group, float('inf'))
        hprio = hdict.get('APpriority', float('inf'))
        if hprio < gprio:
            gprio = hprio

        groups_prio[group] = gprio

    gnames = list(groups_prio.keys())
    gnames.sort(key = lambda g : (groups_prio[g], g) if groups_prio[g] is not None else float('inf'))

    for gname in gnames:
        hset = groups[gname]
        gkw = dict()
        if gname is None:
            apG = ap
        else:
            gdict = groups_kw.get(gname, None)
            if gdict is not None:
                ghelp = gdict.get('about', None)
                if ghelp is not None:
                    gkw['description'] = padding_remove(ghelp)
        apG = ap.add_argument_group(gname, **gkw)
        hlist = list(hset)
        hlist.sort(key = lambda h : (kwdict[h].get("APpriority", float('inf')), h))
        for hname in hlist:
            hdict = kwdict[hname]
            name = hdict.get(hname, hname)

            APshort = hdict.get('APshort', None)
            APlong = '--{}'.format(name)
            APflags = list(hdict.get('APflags', [APlong]))

            if APshort is not None:
                APflags = [APshort] + APflags

            APhide = hdict.get('APhide', False)
            if APhide:
                helptext = argparse.SUPPRESS
            else:
                helptext = hdict.get('about', "<needs doc>").replace('%', '%%')
                helptext = padding_remove(helptext)
            APkw = dict(
                dest = name,
                help = helptext,
            )

            APaction = hdict.get('APaction', None)
            if APaction is not None:
                APkw['action'] = APaction

            APtype = hdict.get('APtype', None)
            if APtype is not None:
                APkw['type'] = APtype

            APrequired = hdict.get('APrequired', None)
            if APrequired is not None:
                APkw['required'] = APrequired

            APnargs = hdict.get('APnargs', None)
            if APnargs is not None:
                APkw['nargs'] = APnargs

            APchoices  = hdict.get('APchoices ', None)
            if APchoices is not None:
                APkw['choices '] = APchoices

            APmetavar = hdict.get('APmetavar', None)
            if APmetavar is not None:
                APkw['metavar'] = APmetavar

            APchoices = hdict.get('APchoices', None)
            if APchoices is not None:
                APkw['choices'] = APchoices

            APconst = hdict.get('APconst', None)
            if APconst is not None:
                APkw['const'] = APconst

            APdefault = hdict.get('APdefault', None)
            if APdefault is not None:
                APkw['default'] = APdefault
            else:
                APkw['default'] = argparse.SUPPRESS

            APpositional = hdict.get('APpositional', False)

            if not APpositional:
                apG.add_argument(*APflags, **APkw)

                #add the aliases to parse, but with their help suppressed
                for aname in hdict.get('aliases', []):
                    APkw['help'] = argparse.SUPPRESS
                    apG.add_argument('--{}'.format(aname), **APkw)
                for aname in hdict.get('aliases_bad', []):
                    APkw['help'] = argparse.SUPPRESS
                    apG.add_argument('--{}'.format(aname), **APkw)
            else:
                dest = APkw.pop('dest')
                apG.add_argument(dest, **APkw)

    return ap
