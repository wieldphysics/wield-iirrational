"""
"""

from wavestate.bunch import Bunch
import numpy as np
import copy
import itertools
import scipy

from . import xfer_algorithms
from . import ss_algorithms


class StateSpaceDense(object):

    def __init__(
        self,
        A        = np.array([[]]),
        B        = np.array([[]]),
        C        = np.array([[]]),
        D        = None,
        E        = None,
        inputs   = None,
        outputs  = None,
        states   = None,
        constr   = None,
        inputsN  = None,
        outputN  = None,
        statesN  = None,
        constrN  = None,
        name     = None,
        N_output = None,
        N_inputs = None,
        N_states = None,
        N_constr = None,
    ):
        self.name = name
        if N_states is None:
            N_states = A.shape[1]
            if N_states == 0:
                N_states = C.shape[1]
            if N_states == 0:
                N_states = B.shape[0]

        if N_constr is None:
            N_constr = A.shape[0]

        if N_inputs is None:
            N_inputs = B.shape[1]
            if N_inputs == 0:
                N_inputs = D.shape[1]

        if N_output is None:
            N_output = C.shape[0]
            if N_output == 0:
                N_output = D.shape[0]

        if E is None:
            if N_constr != N_states:
                raise RuntimeError((
                    "E must be given since the number of "
                    "constraints is different than the number of states"
                ))
            E = np.diag(np.ones(N_states))
        else:
            assert(E.shape == (N_constr, N_states))

        if D is None:
            D = np.zeros((N_output, N_inputs))

        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

        if inputs is None:
            if name is None:
                raise RuntimeError("name must be specified if inputs are not")
            inputs = ["{}.i{}".format(name, idx) for idx in range(N_inputs)]
        if inputsN is None:
            inputsN = np.ones(N_inputs)
        self.inputs = Bunch(
            N        = N_inputs,
            idx2name = inputs,
            idx2N    = np.cumsum([0] + list(inputsN), dtype = int),
            name2idx = {k: i for i, k in enumerate(inputs)},
        )

        if outputs is None:
            if name is None:
                raise RuntimeError("name must be specified if outputs are not")
            outputs = ["{}.o{}".format(name, idx) for idx in range(N_output)]
        if outputN is None:
            outputN = np.ones(N_output)
        self.output = Bunch(
            N        = N_output,
            idx2name = outputs,
            idx2N    = np.cumsum([0] + list(outputN), dtype = int),
            name2idx = {k: i for i, k in enumerate(outputs)},
        )

        if states is None:
            if name is None:
                raise RuntimeError("name must be specified if states are not")
            states = [name]
        if statesN is None:
            statesN = [N_states]
        self.states = Bunch(
            N        = N_states,
            idx2name = states,
            idx2N    = np.cumsum([0] + list(statesN), dtype = int),
            name2idx = {k: i for i, k in enumerate(states)},
        )

        if constr is None:
            if name is None:
                raise RuntimeError("name must be specified if constr are not")
            constr = [name]
        if constrN is None:
            constrN = [N_constr]
        self.constr = Bunch(
            idx2name = constr,
            N        = N_constr,
            idx2N    = np.cumsum([0] + list(constrN), dtype = int),
            name2idx = {k: i for i, k in enumerate(constr)},
        )

        check_SS_group(self.inputs)
        check_SS_group(self.output)
        check_SS_group(self.states)
        check_SS_group(self.constr)
        return

    def group(self, which):
        if which == 'inputs':
            return self.inputs
        elif which == 'output':
            return self.output
        elif which == 'states':
            return self.states
        elif which == 'constr':
            return self.constr
        else:
            raise RuntimeError("State names unrecognized")

    def swap_raw(self, which, rangeA, rangeB):
        A1, A2 = rangeA
        B1, B2 = rangeB
        Sl1 = slice(A1, A2)
        Sl2 = slice(B1, B2)

        def copy1(mat):
            temp = np.copy(mat[Sl1, :])
            mat[Sl1, :] = mat[Sl2, :]
            mat[Sl2, :] = temp

        def copy2(mat):
            temp = np.copy(mat[:, Sl1])
            mat[:, Sl1] = mat[:, Sl2]
            mat[:, Sl2] = temp

        if which == 'inputs':
            copy2(self.B)
            copy2(self.D)
        elif which == 'output':
            copy1(self.C)
            copy1(self.D)
        elif which == 'states':
            copy2(self.C)
            copy2(self.A)
            copy2(self.E)
        elif which == 'constr':
            copy1(self.B)
            copy1(self.A)
            copy1(self.E)
        return

    def percolate_raw(self, which, ranges, keep = True):
        ranges = sorted(ranges)

        if keep:
            def percolate1(mat):
                segments = []
                for r in ranges:
                    A1, A2 = r
                    Sl1 = slice(A1, A2)
                    segments.append(
                        np.copy(mat[Sl1, :])
                    )
                past_r = ranges[0]
                past_idx = past_r[0]
                for r in ranges[1:]:
                    Sl1 = slice(past_r[1], r[0])
                    next_idx = past_idx + r[0] - past_r[1]
                    Sl2 = slice(past_idx, next_idx)
                    mat[Sl2, :] = mat[Sl1, :]
                    past_idx = next_idx
                    past_r = r
                Sl1 = slice(past_r[1], mat.shape[0])
                next_idx = past_idx + mat.shape[0] - past_r[1]
                Sl2 = slice(past_idx, next_idx)
                mat[Sl2, :] = mat[Sl1, :]
                for seg in segments:
                    next_idx = past_idx + seg.shape[0]
                    Sl2 = slice(past_idx, next_idx)
                    mat[Sl2, :] = seg
                    past_idx = next_idx

            def percolate2(mat):
                segments = []
                for r in ranges:
                    A1, A2 = r
                    Sl1 = slice(A1, A2)
                    segments.append(
                        np.copy(mat[:, Sl1])
                    )
                past_r = ranges[0]
                past_idx = past_r[0]
                for r in ranges[1:]:
                    Sl1 = slice(past_r[1], r[0])
                    next_idx = past_idx + r[0] - past_r[1]
                    Sl2 = slice(past_idx, next_idx)
                    mat[:, Sl2] = mat[:, Sl1]
                    past_idx = next_idx
                    past_r = r
                Sl1 = slice(past_r[1], mat.shape[1])
                next_idx = past_idx + mat.shape[1] - past_r[1]
                Sl2 = slice(past_idx, next_idx)
                mat[:, Sl2] = mat[:, Sl1]
                for seg in segments:
                    next_idx = past_idx + seg.shape[1]
                    Sl2 = slice(past_idx, next_idx)
                    mat[:, Sl2] = seg
                    past_idx = next_idx
        else:
            def percolate1(mat):
                past_r = ranges[0]
                past_idx = past_r[0]
                for r in ranges[1:]:
                    Sl1 = slice(past_r[1], r[0])
                    next_idx = past_idx + r[0] - past_r[1]
                    Sl2 = slice(past_idx, next_idx)
                    mat[Sl2, :] = mat[Sl1, :]
                    past_idx = next_idx
                    past_r = r
                Sl1 = slice(past_r[1], mat.shape[0])
                next_idx = past_idx + mat.shape[0] - past_r[1]
                Sl2 = slice(past_idx, next_idx)
                mat[Sl2, :] = mat[Sl1, :]

            def percolate2(mat):
                past_r = ranges[0]
                past_idx = past_r[0]
                for r in ranges[1:]:
                    Sl1 = slice(past_r[1], r[0])
                    next_idx = past_idx + r[0] - past_r[1]
                    Sl2 = slice(past_idx, next_idx)
                    mat[:, Sl2] = mat[:, Sl1]
                    past_idx = next_idx
                    past_r = r
                Sl1 = slice(past_r[1], mat.shape[1])
                next_idx = past_idx + mat.shape[1] - past_r[1]
                Sl2 = slice(past_idx, next_idx)
                mat[:, Sl2] = mat[:, Sl1]

        if which == 'inputs':
            percolate2(self.B)
            percolate2(self.D)
        elif which == 'output':
            percolate1(self.C)
            percolate1(self.D)
        elif which == 'states':
            percolate2(self.C)
            percolate2(self.A)
            percolate2(self.E)
        elif which == 'constr':
            percolate1(self.B)
            percolate1(self.A)
            percolate1(self.E)
        else:
            raise RuntimeError("Unrecognized 'which' argument")
        return

    def names_collect(
            self,
            which,
            to,
            fr = None,
    ):
        """
        Rename states

        which is the type, must be one of "inputs, output, states, constr"
        """
        G = self.group(which)

        if fr is None:
            #completely rename
            G.idx2name.clear()
            G.idx2name = [to]
            G.idx2N = np.asarray([0, G.N])
            G.name2idx = {to : 0}
        else:
            raise NotImplementedError()

    def names_split(
            self,
            which,
            to,
            fr,
    ):
        """
        split state names

        which is the type, must be one of "inputs, output, states, constr"
        """
        #G = self.group(which)
        raise NotImplementedError()

    def names_change(
            self,
            which,
            to,
            fr,
    ):
        """
        Rename states

        which is the type, must be one of "inputs, output, states, constr"
        """
        G = self.group(which)

        idx = G.name2idx.pop(fr)
        G.name2idx[to] = idx
        G.idx2name[idx] = to
        return

    def names_get(
        self,
        which,
        name,
    ):
        """
        Get the index and size for name
        """
        G = self.group(which)

        idx = G.name2idx[name]
        idxN = G.idx2N[idx]
        idxN2 = G.idx2N[idx + 1]
        return idxN, idxN2 - idxN

    def __getitem__(self, key):
        which, name = key
        return self.names_get(which, name)

    def outputs_delete(self, names):
        ranges = []
        Nlist = list(self.output.idx2N[1:] - self.output.idx2N[:-1])
        idx_del = []
        for name in names:
            idx = self.output.name2idx.pop(name)
            self.output.idx2N[idx]
            ranges.append(
                (self.output.idx2N[idx], self.output.idx2N[idx + 1])
            )
            idx_del.append(idx)

        idx_del = sorted(idx_del)
        #must go highest to lowest due to index deletion
        for idx in reversed(idx_del):
            del Nlist[idx]
            del self.output.idx2name[idx]
        #now need to fix the name2idx list since the idx have changed
        #reuse the last idx of the last iteration
        for idx in range(idx, len(self.output.idx2name)):
            name = self.output.idx2name[idx]
            self.output.name2idx[name] = idx

        self.output.idx2N = np.cumsum([0] + list(Nlist), dtype = int)
        self.output.N = self.output.idx2N[-1]

        self.percolate_raw('output', ranges, keep = False)
        Ntotal = 0
        for r in ranges:
            Ntotal += r[1] - r[0]

        self.D = self.D[:-Ntotal, :]
        self.C = self.C[:-Ntotal, :]
        return

    def iter(
        self,
        which,
    ):
        """
        Rename states

        which is the type, must be one of "inputs, output, states, constr"
        """
        yield

    def in2out(self):
        """
        Modify the statespace so that the inputs become outputs for
        output binding. Adds one space per input
        """
        N_states = self.states.N + self.inputs.N
        N_output = self.output.N + self.inputs.N

        Zb = np.zeros_like(self.B)
        Oi = np.diag(np.ones(self.inputs.N))

        A = np.block([
            [self.A, self.B],
        ])

        E = np.block([
            [self.E, Zb],
        ])

        C = np.block([
            [self.C, self.D],
            [Zb.T, Oi]
        ])

        D = np.array([]).reshape((N_output, 0))
        B = np.array([]).reshape((N_states, 0))

        new = self.__class__.__new__(self.__class__)
        new.A = A
        new.E = E
        new.B = B
        new.C = C
        new.D = D

        new.constr = copy.deepcopy(self.constr)
        new.inputs = Bunch(
            idx2name = [],
            idx2N    = np.asarray([0]),
            name2idx = {},
            N        = 0,
        )
        outlen = len(self.output.idx2name)
        outN   = self.output.N
        name2idx = dict(self.output.name2idx)
        for iname, idx in self.inputs.name2idx.items():
            assert(iname not in name2idx)
            name2idx[iname] = idx + outlen
        new.output = Bunch(
            idx2name = self.output.idx2name + self.inputs.idx2name,
            idx2N    = np.concatenate([self.output.idx2N, outN + self.inputs.idx2N[1:]]),
            name2idx = name2idx,
            N        = self.output.N + self.inputs.N,
        )

        stlen = len(self.states.idx2name)
        stN   = self.states.N
        name2idx = dict(self.states.name2idx)
        for iname, idx in self.inputs.name2idx.items():
            assert(iname not in name2idx)
            name2idx[iname] = idx + stlen
        new.states = Bunch(
            idx2name = self.states.idx2name + self.inputs.idx2name,
            idx2N    = np.concatenate([self.states.idx2N, stN + self.inputs.idx2N[1:]]),
            name2idx = name2idx,
            N        = self.states.N + self.inputs.N,
        )
        check_SS_group(new.inputs)
        check_SS_group(new.output)
        check_SS_group(new.states)
        check_SS_group(new.constr)
        new.name = self.name
        return new

    def xfer(self, F_Hz, iname, oname, isub = 0, osub = 0):
        idx_i = self.inputs.name2idx[iname]
        idx_o = self.output.name2idx[oname]
        return xfer_algorithms.ss2xfer(
            A       = self.A,
            B       = self.B,
            C       = self.C,
            D       = self.D,
            E       = self.E,
            F_Hz    = F_Hz,
            idx_in  = idx_i + isub,
            idx_out = idx_o + osub,
        )

    def constrain(self, name, F, outputs):
        """
        Modifies DSS to have an additional constraint given

        TODO, add ability to have many constraint names
        """
        return self.constraints(
            [
                name,
                (F, outputs)
            ]
        )

    def constraints(self, name, bond_list):
        """
        Modifies DSS to have an additional constraint given

        TODO, add ability to have many constraint names
        """
        Nplus = 0
        all_outs = set()
        Aadd = []
        Dadd = []

        def update(F, outputs):
            """
            Internal function to generate constraint equations
            """
            nonlocal Nplus

            F = np.asarray(F)
            if len(F.shape) == 1:
                F = F.reshape(1, -1)

            Cmat = []
            Dmat = []
            for oname in outputs:
                oidx = self.output.name2idx[oname]
                oN1 = self.output.idx2N[oidx]
                oN2 = self.output.idx2N[oidx + 1]
                Cmat.append(self.C[oN1:oN2, :])
                Dmat.append(self.D[oN1:oN2, :])

            Cmat = np.vstack(Cmat)
            Dmat = np.vstack(Dmat)
            Aconstr = F @ Cmat
            Dconstr = F @ Dmat
            Aadd.append(Aconstr)
            Dadd.append(Dconstr)
            Nplus += F.shape[0]

        def commit(name):
            """
            Internal function to commit constraint equations to the statespace
            """
            nonlocal Aadd, Dadd
            if len(Aadd) is None:
                return
            Aadd = np.vstack(Aadd)
            Dadd = np.vstack(Dadd)

            assert(name not in self.constr.name2idx)

            A = np.block([
                [self.A,  ],
                [Aadd, ],
            ])

            E = np.block([
                [self.E,                 ],
                [np.zeros_like(Aadd), ],
            ])

            B = np.block([
                [self.B,  ],
                [Dadd, ]
            ])

            self.constr.name2idx[name] = len(self.constr.idx2name)
            self.constr.idx2name.append(name)
            self.constr.idx2N = np.concatenate([self.constr.idx2N, [self.constr.idx2N[-1] + Nplus]])
            self.constr.N += Nplus
            self.A = A
            self.E = E
            self.B = B
            assert(self.A.shape[0] == self.constr.N)
            assert(self.E.shape[0] == self.constr.N)
            check_SS_group(self.constr)

        for bond in bond_list:
            if isinstance(bond, str):
                commit(name)
                name = bond
                continue

            F = bond[0]

            if isinstance(F, str):
                F = F.lower()
                if F == 'bind':
                    outputs_grp = bond[1:]
                    for outputs in outputs_grp:
                        all_outs.update(outputs)
                        constr = np.array([[-1, 1]])
                        for pair in itertools.combinations(outputs, 2):
                            update(constr, pair)
                elif F == 'flow':
                    outputs_grp = bond[1:]
                    for outputs in outputs_grp:
                        all_outs.update(outputs)
                        update(np.ones(len(outputs)), outputs)
                elif F == 'zero':
                    outputs_grp = bond[1:]
                    for outputs in outputs_grp:
                        all_outs.update(outputs)
                        for output in outputs:
                            update([[1]], outputs)
                elif F == 'sum_into':
                    F, output_to, outputs_fr = bond
                    all_outs.add(output_to)
                    all_outs.update(outputs_fr)
                    constr = np.ones(len(outputs_fr) + 1)
                    constr[0] = -1
                    update(constr, (output_to,) + tuple(outputs_fr))
                else:
                    raise RuntimeError("Unrecognized bond type")
            else:
                #then F is a matrix to be used with each set of outputs in the group
                F, outputs_grp = bond
                for outputs in outputs_grp:
                    all_outs.update(outputs)
                    update(F, outputs)
        #now do a final commit since there can be no more names
        commit(name)

        return all_outs

    @classmethod
    def join(cls, name, SSs):
        self = cls.__new__(cls)
        self.name = name

        self.inputs = Bunch(
            idx2name = [],
            idx2N    = [0],
            name2idx = {},
        )
        self.output = Bunch(
            idx2name = [],
            idx2N    = [0],
            name2idx = {},
        )
        self.states = Bunch(
            idx2name = [],
            idx2N    = [0],
            name2idx = {},
        )
        self.constr = Bunch(
            idx2name = [],
            idx2N    = [0],
            name2idx = {},
        )

        ss_seq = []
        for ss in SSs:
            seq_inputs = []
            seq_output = []
            seq_states = []
            seq_constr = []
            Gtransfer(ss.inputs, self.inputs, seq_inputs)
            Gtransfer(ss.output, self.output, seq_output)
            Gtransfer(ss.states, self.states, seq_states)
            Gtransfer(ss.constr, self.constr, seq_constr)
            #print((ss, seq_inputs, seq_output, seq_states, seq_constr))
            ss_seq.append(
                (ss, seq_inputs, seq_output, seq_states, seq_constr)
            )

        self.inputs.idx2N = np.asarray(self.inputs.idx2N)
        self.output.idx2N = np.asarray(self.output.idx2N)
        self.states.idx2N = np.asarray(self.states.idx2N)
        self.constr.idx2N = np.asarray(self.constr.idx2N)
        self.inputs.N = self.inputs.idx2N[-1]
        self.output.N = self.output.idx2N[-1]
        self.states.N = self.states.idx2N[-1]
        self.constr.N = self.constr.idx2N[-1]

        check_SS_group(self.inputs)
        check_SS_group(self.output)
        check_SS_group(self.states)
        check_SS_group(self.constr)

        A = np.zeros((self.constr.N, self.states.N))
        E = np.zeros((self.constr.N, self.states.N))
        B = np.zeros((self.constr.N, self.inputs.N))
        C = np.zeros((self.output.N, self.states.N))
        D = np.zeros((self.output.N, self.inputs.N))

        for (ss, seq_inputs, seq_output, seq_states, seq_constr) in ss_seq:
            for (Nc1f, Nc2f, Nc1t, Nc2t) in seq_constr:
                for (Ns1f, Ns2f, Ns1t, Ns2t) in seq_states:
                    A[Nc1t : Nc2t, Ns1t : Ns2t] = ss.A[Nc1f : Nc2f, Ns1f : Ns2f]
                    E[Nc1t : Nc2t, Ns1t : Ns2t] = ss.E[Nc1f : Nc2f, Ns1f : Ns2f]

                for (Ni1f, Ni2f, Ni1t, Ni2t) in seq_inputs:
                    B[Nc1t : Nc2t, Ni1t : Ni2t] = ss.B[Nc1f : Nc2f, Ni1f : Ni2f]

            for (No1f, No2f, No1t, No2t) in seq_output:
                for (Ns1f, Ns2f, Ns1t, Ns2t) in seq_states:
                    C[No1t : No2t, Ns1t : Ns2t] = ss.C[No1f : No2f, Ns1f : Ns2f]

                for (Ni1f, Ni2f, Ni1t, Ni2t) in seq_inputs:
                    D[No1t : No2t, Ni1t : Ni2t] = ss.D[No1f : No2f, Ni1f : Ni2f]

        self.A = A
        self.E = E
        self.B = B
        self.C = C
        self.D = D

        check_SS_group(self.inputs)
        check_SS_group(self.output)
        check_SS_group(self.states)
        check_SS_group(self.constr)
        return self

    @classmethod
    def chain(cls, name, SSs):
        self = cls.__new__(cls)
        self.name = name

        self.inputs = copy.deepcopy(SSs[0].inputs)
        self.output = copy.deepcopy(SSs[-1].output)

        self.states = Bunch(
            idx2name = [],
            idx2N    = [0],
            name2idx = {},
        )
        self.constr = Bunch(
            idx2name = [],
            idx2N    = [0],
            name2idx = {},
        )

        ss_seq = []
        for ss in SSs:
            seq_states = []
            seq_constr = []
            Gtransfer(ss.states, self.states, seq_states)
            Gtransfer(ss.constr, self.constr, seq_constr)
            #print((ss, seq_states, seq_constr))
            ss_seq.append(
                (ss, seq_states, seq_constr)
            )

        self.states.idx2N = np.asarray(self.states.idx2N)
        self.constr.idx2N = np.asarray(self.constr.idx2N)
        self.states.N = self.states.idx2N[-1]
        self.constr.N = self.constr.idx2N[-1]

        check_SS_group(self.states)
        check_SS_group(self.constr)

        A = np.zeros((self.constr.N, self.states.N))
        E = np.zeros((self.constr.N, self.states.N))
        B = np.zeros((self.constr.N, self.inputs.N))
        C = np.zeros((self.output.N, self.states.N))

        D = None
        for idx_ss, (ss, seq_states, seq_constr) in enumerate(ss_seq):
            for (Nc1f, Nc2f, Nc1t, Nc2t) in seq_constr:
                for (Ns1f, Ns2f, Ns1t, Ns2t) in seq_states:
                    A[Nc1t : Nc2t, Ns1t : Ns2t] = ss.A[Nc1f : Nc2f, Ns1f : Ns2f]
                    E[Nc1t : Nc2t, Ns1t : Ns2t] = ss.E[Nc1f : Nc2f, Ns1f : Ns2f]
            if idx_ss > 0:
                B_ud = ss.B
                idx_down = idx_ss - 1
                while True:
                    (ss_down, seq_states_down, seq_constr_down) = ss_seq[idx_down]
                    A_ud = B_ud @ ss_down.C
                    for (Ns1f, Ns2f, Ns1t, Ns2t) in seq_states_down:
                        for (Nc1f, Nc2f, Nc1t, Nc2t) in seq_constr:
                            A[Nc1t : Nc2t, Ns1t : Ns2t] = A_ud[Nc1f : Nc2f, Ns1f : Ns2f]
                    if idx_down == 0:
                        break
                    B_ud = B_ud @ ss_down.D
                    idx_down -= 1

                B_into = ss.B @ D
                for (Nc1f, Nc2f, Nc1t, Nc2t) in seq_constr:
                    B[Nc1t : Nc2t, :] = B_into[Nc1f : Nc2f, :]
                D = ss.D @ D
            else:
                for (Nc1f, Nc2f, Nc1t, Nc2t) in seq_constr:
                    B[Nc1t : Nc2t, :] = ss.B[Nc1f : Nc2f, :]
                D = ss.D

        (ss, seq_states, seq_constr) = ss_seq[-1]
        for (Ns1f, Ns2f, Ns1t, Ns2t) in seq_states:
            C[:, Ns1t : Ns2t] = ss.C[:, Ns1f : Ns2f]
        D_rev = ss.D
        #now loop down through the D matrices
        for (ss, seq_states, seq_constr) in ss_seq[-2::-1]:
            C_into = D_rev @ ss.C
            for (Ns1f, Ns2f, Ns1t, Ns2t) in seq_states:
                C[:, Ns1t : Ns2t] = C_into[:, Ns1f : Ns2f]
            D_rev = D_rev @ ss.D

        self.A = A
        self.E = E
        self.B = B
        self.C = C
        self.D = D

        check_SS_group(self.states)
        check_SS_group(self.constr)
        return self

    def inverse(
            self,
            name = None,
            augment_states = 'output',
            augment_constr = 'new',
            method = 'DSS',
            name_transformer = lambda n : n
    ):
        new = self.__class__.__new__(self.__class__)

        if self.output.N != self.inputs.N:
            raise RuntimeError("Cannot invert")
        N = self.output.N

        new.inputs = copy.deepcopy(self.output)
        new.output = copy.deepcopy(self.inputs)
        new.states = copy.deepcopy(self.states)
        new.constr = copy.deepcopy(self.constr)

        if name is None:
            name = self.name
        new.name = name

        def inject(Gfr, Gto):
            for idx, name in enumerate(Gfr.idx2name):
                name_alt = name_transformer(name)
                idx_new = len(Gto.idx2name)
                idx_check = Gto.name2idx.setdefault(name_alt, idx_new)
                if idx_check != idx_new:
                    raise RuntimeError("Name conflict on inversion")
                Gto.idx2name.append(name_alt)
            Gto.idx2N = np.concatenate([Gto.idx2N, Gto.idx2N[-1] + Gfr.idx2N[1:]])
            Gto.N = Gto.N + Gfr.N

        def create(name, N, Gto):
            idx_new = len(Gto.idx2name)
            idx_check = Gto.name2idx.setdefault(name, idx_new)
            if idx_check != idx_new:
                raise RuntimeError("Name conflict on inversion")
            Gto.idx2name.append(name)
            Gto.idx2N = np.concatenate([Gto.idx2N, [Gto.idx2N[-1] + N]])
            Gto.N = Gto.N + N

        if augment_states == 'inputs':
            inject(new.inputs, new.states)
        elif augment_states == 'output':
            inject(new.output, new.states)
        elif augment_states == 'new':
            create('{}:inv'.format(self.name), N, new.states)
        else:
            raise RuntimeError("augment_names must be one of 'inputs', 'output', 'new'")

        if augment_constr == 'inputs':
            inject(new.inputs, new.constr)
        elif augment_constr == 'output':
            inject(new.output, new.constr)
        elif augment_constr == 'new':
            create('{}:inv'.format(self.name), N, new.constr)
        else:
            raise RuntimeError("augment_names must be one of 'inputs', 'output', 'new'")

        A = np.zeros((new.constr.N, new.states.N))
        E = np.zeros((new.constr.N, new.states.N))
        B = np.zeros((new.constr.N, new.inputs.N))
        C = np.zeros((new.output.N, new.states.N))
        D = np.zeros((new.output.N, new.inputs.N))

        A[:self.constr.N, :self.states.N] = self.A
        A[:self.constr.N, self.states.N:] = self.B
        A[self.constr.N:, :self.states.N] = self.C
        A[self.constr.N:, self.states.N:] = self.D
        E[:self.constr.N, :self.states.N] = self.E
        B[self.constr.N:, :self.inputs.N] = -1
        C[:self.output.N, self.states.N:] = 1

        new.A = A
        new.E = E
        new.B = B
        new.C = C
        new.D = D

        return new

    def reduce(self, constr = None, states = None, tol = 1e-10):
        #from icecream import ic

        #properly dealing with names is going to be painful
        if states is None:
            states_nZranges = [(0, self.states.N)]
            self.names_collect('states', to = 'reduced')
        else:
            raise RuntimeError("Does not yet support specifying states for reduction")

        if constr is None:
            constr_nZranges = [(0, self.constr.N)]
            self.names_collect('constr', to = 'reduced')
        else:
            raise RuntimeError("Does not yet support specifying constraints for reduction")

        statemask = None
        constr_Zs = []
        for constrR in constr_nZranges:
            constrSl = slice(*constrR)
            for statesR in states_nZranges:
                statesSl = slice(*statesR)
                #TODO, make the zero check a touch smarter for symbolics
                mask = ~np.any(self.E[constrSl, statesSl], axis = 1)
                #ic(mask)
                if statemask is None:
                    statemask = mask
                else:
                    statemask = statemask & mask
            constr_Zs.append(constrR[0] + np.argwhere(statemask).reshape(-1))
        constr_Zranges = argwhere2ranges(np.concatenate(constr_Zs))
        constr_nZranges = anti_ranges(self.constr.N, constr_Zranges)

        statemask = None
        states_Zs = []
        for statesR in states_nZranges:
            statesSl = slice(*statesR)
            for constrR in constr_nZranges:
                constrSl = slice(*constrR)
                #TODO, make the zero check a touch smarter for symbolics
                mask = ~np.any(self.E[constrSl, statesSl], axis = 0)
                #ic(mask)
                if statemask is None:
                    statemask = mask
                else:
                    statemask = statemask & mask
            states_Zs.append(statesR[0] + np.argwhere(statemask).reshape(-1))
        states_Zranges = argwhere2ranges(np.concatenate(states_Zs))
        states_nZranges = anti_ranges(self.states.N, states_Zranges)

        for R in constr_Zranges:
            assert(not np.any(self.E[slice(*R), :]))
        for R in states_Zranges:
            assert(not np.any(self.E[:, slice(*R)]))

        #now rearrange
        Nst = sum(R1 - R0 for R0, R1 in states_Zranges)
        Nco = sum(R1 - R0 for R0, R1 in constr_Zranges)
        Ist = self.states.N - Nst
        Ico = self.constr.N - Nco
        self.percolate_raw('states', states_Zranges)
        self.percolate_raw('constr', constr_Zranges)
        #print(Nst, Nco)
        assert(Nst == Nco)
        assert(not np.any(self.E[Ico:, :]))
        assert(not np.any(self.E[:, Ist:]))

        A22 = self.A[Ico:, Ist:]
        U, S, V = scipy.linalg.svd(A22)

        #now reorder for smallest first
        S = S[::-1]
        Nlast = len(S) - np.searchsorted(S, S[-1] * tol)
        #print("NLAST", Nlast)
        #redefines constr
        U = U[:, ::-1]
        Ui = U.T
        #redefines states
        V = V[::-1, :]
        Vi = V.T

        #ic(S)
        #ic(np.diagonal(mat))
        self.B[Ico:, :] = Ui @ self.B[Ico:, :]
        self.C[:, Ist:] = self.C[:, Ist:] @ Vi
        #A22
        #one might expect to use diag(S) here, but
        #the svd is still unstable, so we need to use only orthogonal
        #transformations
        if Nlast == len(S):
            self.A[Ico:, Ist:] = np.diag(S)
        else:
            self.A[Ico:, Ist:] = Ui @ self.A[Ico:, Ist:] @ Vi
        #A12
        self.A[:Ico, Ist:] = self.A[:Ico, Ist:] @ Vi
        #A21
        self.A[Ico:, :Ist] = Ui @ self.A[Ico:, :Ist]

        #now, select the last few rows that are OK
        #and reduce them down
        S = S[-Nlast:]
        Ist = self.states.N - Nlast
        Ico = self.constr.N - Nlast

        A22 = self.A[Ico:, Ist:]
        if False:
            #this secondary SVD is useful to improve
            #numerical stability yet further, but appears to be unnecessary
            U, S, V = scipy.linalg.svd(A22)
            #now reorder for smallest first
            #redefines constr
            Ui = U.T
            #redefines states
            Vi = V.T
            #ic(S)
            self.B[Ico:, :] = Ui @ self.B[Ico:, :]
            self.C[:, Ist:] = self.C[:, Ist:] @ Vi

            self.B[Ico:, :] = Ui @ self.B[Ico:, :]
            self.C[:, Ist:] = self.C[:, Ist:] @ Vi
            #A12
            self.A[:Ico, Ist:] = self.A[:Ico, Ist:] @ Vi
            #A21
            self.A[Ico:, :Ist] = Ui @ self.A[Ico:, :Ist]

        A12mod = self.A[:Ico, Ist:] / S.reshape(1, -1)
        self.A[:Ico, :Ist] -= A12mod @ self.A[Ico:, :Ist]
        self.B[:Ico, :] -= A12mod @ self.B[Ico:, :]

        C2mod = self.C[:, Ist:] / S.reshape(1, -1)
        self.C[:, :Ist] -= C2mod @ self.A[Ico:, :Ist]
        self.D -= C2mod @ self.B[Ico:, :]

        self.C = self.C[:, :Ist]
        self.B = self.B[:Ico, :]
        self.A = self.A[:Ico, :Ist]
        self.E = self.E[:Ico, :Ist]

        #works because of the state name reduction above
        #TODO, do the proper splitting/absorption of the reduced states
        self.constr.N -= Nlast
        self.states.N -= Nlast
        self.constr.idx2N[-1] -= Nlast
        self.states.idx2N[-1] -= Nlast

        return

    def controllable_staircase(self, constr = None, states = None, tol = 1e-10):
        #properly dealing with names is going to be painful
        if states is None:
            states_nZranges = [(0, self.states.N)]
            self.names_collect('states', to = 'reduced')
        else:
            raise RuntimeError("Does not yet support specifying states for reduction")

        if constr is None:
            constr_nZranges = [(0, self.constr.N)]
            self.names_collect('constr', to = 'reduced')
        else:
            raise RuntimeError("Does not yet support specifying constraints for reduction")

        A, B, C, D, E = ss_algorithms.controllable_staircase(
            self.A,
            self.B,
            self.C,
            self.D,
            self.E,
            tol = tol
        )
        self.A = A
        self.E = E
        self.B = B
        self.C = C
        self.D = D
        return

def check_SS_group(G):
    assert(len(G.idx2N) == 1 + len(G.idx2name))
    assert(G.N == G.idx2N[-1])
    assert(G.idx2N is np.asarray(G.idx2N))
    assert(isinstance(G.idx2name, list))


def Gtransfer(Gfr, Gto, seq):
    idx_fr_st = None
    idx_to_st = None
    idx_fr_end = None
    idx_to_end = None

    def do_seq():
        nonlocal idx_fr_st, idx_to_st, idx_fr_end, idx_to_end
        if idx_to_st is None:
            return
        #if idx_fr_st != idx_fr_end and idx_to_st != idx_to_end:
        N_fr_st = Gfr.idx2N[idx_fr_st]
        N_fr_end = Gfr.idx2N[idx_fr_end + 1]

        N_to_st = Gto.idx2N[idx_to_st]
        N_to_end = Gto.idx2N[idx_to_end + 1]
        assert(N_to_end - N_to_st == N_fr_end - N_fr_st)
        seq.append((N_fr_st, N_fr_end, N_to_st, N_to_end))

    for idx_fr, iname in enumerate(Gfr.idx2name):
        if idx_fr_st is None:
            idx_fr_st = idx_fr_end = idx_fr
        idx_to = Gto.name2idx.get(iname, None)
        N = Gfr.idx2N[idx_fr + 1] - Gfr.idx2N[idx_fr]
        if idx_to is None:
            idx_to = len(Gto.idx2name)
            Gto.idx2name.append(iname)
            Gto.idx2N.append(Gto.idx2N[-1] + N)
            Gto.name2idx[iname] = idx_to
            idx_fr_end = idx_fr
            if idx_to_st is None:
                idx_to_st = idx_to
            idx_to_end = idx_to
        else:
            #check size
            N2 = Gto.idx2N[idx_to + 1] - Gto.idx2N[idx_to]
            assert(N == N2)
            do_seq()
            idx_fr_st = idx_fr_end = idx_fr
            idx_to_st = idx_to_end = idx_to
    do_seq()


def argwhere2ranges(wheres):
    delta = (wheres[1:] - wheres[:-1] != 1)
    ranges = []
    idx_p = wheres[0]
    for idx in np.argwhere(delta).reshape(-1):
        idx_n = wheres[idx] + 1
        ranges.append((idx_p, idx_n))
        idx_p = wheres[idx + 1]
    ranges.append((idx_p, wheres[-1] + 1))
    return ranges


def anti_ranges(N, ranges):
    idx_p = 0
    anti_ranges = []
    for R in ranges:
        idx_n = R[0]
        if idx_p != idx_n:
            anti_ranges.append((idx_p, idx_n))
        idx_p = R[1]
    idx_n = N
    if idx_p != idx_n:
        anti_ranges.append((idx_p, idx_n))
    return anti_ranges

