"""
CABSDock module for contact map analysis of trajectories.

Created on 4 June 2017 by Tymoteusz hert Oleniecki.
"""
import numpy
import operator

from cabsDock.plots import mk_histo
from cabsDock.utils import _chunk_lst

import matplotlib.pyplot
import matplotlib.ticker
import matplotlib.colors


class ContactMapFactory(object):
    def __init__(self, chains1, chains2, temp):
        """Builder for ContactMap.

        Arguments:
        chains1 -- list or str; chars for 1st chain(s).
        chains2 -- list or str; chars for 2nd chain(s).
        temp -- cabsDock.atom.Atoms instance containing both given chains.

        """
        chs = {}
        for n, i in enumerate(temp.atoms):
            chs.setdefault(i.chid, []).append(n)
        self.dims = (sum(map(len, [chs.get(ch1, []) for ch1 in chains1])), sum(map(len, [chs.get(ch2, []) for ch2 in chains2])))
        self.inds1 = reduce(operator.add, [chs.get(i, []) for i in chains1])
        self.inds2 = reduce(operator.add, [chs.get(i, []) for i in chains2])
        self.ats1 = [temp.atoms[i] for i in self.inds1]
        self.ats2 = [temp.atoms[i] for i in self.inds2]

    def mk_cmap(self, traj, thr, frames=None, replicas=None):
        """Creates map of contacts between two given chains.

        Arguments:
        traj -- numpy.array of propper shape, i.e. Nreplicas x Nframes x Natoms x 3.
        thr -- float; threshold for side chain distance contact.
        frames -- tuple of ints; indexes of frames to be taken. All frames are taken by default.
        replicas -- list of replicas' indexes to be taken. All replicas are taken by default.

        Returns list of ContactMap for each replica in trajectory.
        """
        
        if replicas is None:
            replicas = xrange(traj.shape[0])
        if frames is None:
            fstf = 0
            frames = slice(1, None)
        else:
            frames = list(frames)
            fstf = frames.pop(0)
        resl = []
        for rep in traj[replicas,]:
            cmtx = self.mk_cmtx(self.mk_dmtx(rep[fstf]), thr)
            nframes = 1
            for fra in rep[frames,]:
                ncmtx = self.mk_cmtx(self.mk_dmtx(fra), thr)
                cmtx += ncmtx
                nframes += 1
            resl.append(ContactMap(cmtx, [i.fmt() for i in self.ats1], [i.fmt() for i in self.ats2], nframes))
        return resl

    def mk_cmtx(self, mtx, thr):
        """Returns boolean numpy.array of contacts from given mtx of distances.

        Arguments:
        mtx -- numpy.array of distances between atoms.
        thr -- thresohld below which atoms are in contact.
        """
        return numpy.clip(numpy.sign(-mtx + thr), 0, 1)

    def mk_dmtx(self, vec):
        """Returns 2D numpy.array of distances between atoms from vector of coords.

        Arguments:
        vec -- slice of trajectory.

        Calculates distances between atoms from given chains.
        """
        m1 = vec[self.inds1,].reshape(-1, 1, 3)[:, self.dims[1] * (0,)]
        m2 = vec[self.inds2,].reshape(1, -1, 3)[self.dims[0] * (0,), :]
        mtx = m1 - m2
        return (mtx * mtx).sum(axis=2) ** .5


class ContactMap(object):
    def __init__(self, mtx, nms1, nms2, n):
        """Contact map init.

        Arguments:
        mtx -- 2D numpy.array of distances between (pseudo)atoms.
        atoms1, atoms2 -- cabsDock.atom.Atoms instance; template for cmap.
        n -- number of frames.
        """
        self.cmtx = mtx
        self.s1 = nms1
        self.s2 = nms2
        self.n = n

    def zero_diagonal(self):
        numpy.fill_diagonal(self.cmtx, 0)

    def save_fig(self, fname, fmt='svg', norm_n=False, break_long_x=50):
        """Saves cmap as matrix plot.

        Arguments:
        fname -- str; file name.
        fmt -- str; by default 'svg'. See matplotlib.pyplot.savefig for more inforamtion.
        norm_n -- bool; if True cmap will be normalized by cmap number of frames.
        break_long_x -- int; 50 by default. If not set to 0 -- chunks long x axis into fragments of given length.
        """
        #TODO: taking min and max out of data + [50] should be done once
        wdh_cnst = 50
        wdth = self.cmtx.shape[0]
        chunks = _chunk_lst(range(wdth), break_long_x) if break_long_x else [range(wdth)]
        lngst = max(map(len, chunks))
        label_size = min(lngst, wdh_cnst) / float(wdh_cnst) * 10
        size = (min(lngst, wdh_cnst) / 5., len(chunks) * min(len(self.s2), wdh_cnst) / 5. + 2)
        fig = matplotlib.pyplot.figure(figsize=size)
        grid_wdth = min(lngst, wdh_cnst)
        grid = matplotlib.pyplot.GridSpec(len(chunks) + 1, 1, height_ratios=([len(self.s2) for i in chunks] + [len(self.s2) * .25]))
        vmax = self.n if norm_n else numpy.max(self.cmtx)
        if vmax < 5:
            vmax = 1 if norm_n else 5
        colors = matplotlib.colors.LinearSegmentedColormap.from_list('bambi',
                zip([0., .01, .1, .4, .7, 1.], ['#ffffff', '#f2d600', '#e80915', '#666666', '#4b8f24', '#000000']))

        for n, chunk in enumerate(chunks):
            sfig = matplotlib.pyplot.subplot(grid[n : n + 1, 0])
            mtx = self.cmtx.T[:,chunk]
            if mtx.shape[1] < grid_wdth:
                zrs = numpy.zeros((mtx.shape[0], grid_wdth))
                zrs[:, :len(chunk)] = mtx
                mtx = zrs
            sfig.matshow(
                mtx,
                cmap=colors,
                vmin=0.,
                vmax=vmax,
                )

            if sfig.get_data_ratio() < 1.:
                aratio = float.__div__(*map(float, mtx.shape))
                sfig.set_aspect(sfig.get_data_ratio() ** -1 * aratio)
            settings = (
                        (list(numpy.array(self.s1)[chunk,]), len(chunk), sfig.set_xticks, sfig.set_xticklabels),
                        (self.s2, len(self.s2), sfig.set_yticks, sfig.set_yticklabels),
                        )
            for lbls, n_tcks, tck_loc, tck_lab in settings:
                nloc = break_long_x if break_long_x else wdh_cnst
                ntcks = min(nloc, n_tcks)
                inds = numpy.linspace(0, len(lbls) -1, ntcks).astype(int)
                locator = matplotlib.ticker.MultipleLocator(len(lbls) / ntcks)
                tck_loc(inds)
                tck_lab(list(numpy.array(lbls)[inds,]))

            sfig.tick_params(labelsize=label_size, bottom=True, top=False, labelbottom=True,labeltop=False)
            for tick in sfig.get_xticklabels():
                tick.set_rotation(90)

        if norm_n:
            n_ticks = 5
            vls = numpy.linspace(0., 1., 256)
            tcks = numpy.linspace(0., 1., n_ticks)
            vmax = 1.
        else:
            vls = numpy.arange(0, int(vmax), 1)
            n_ticks = 15 if len(vls) > 30 else int(vmax)
            tcks = numpy.linspace(0, len(vls) - 1, n_ticks).astype(int)
        tcks_loc = numpy.linspace(0., len(vls) - 1, n_ticks).astype(int)
        ax2 = matplotlib.pyplot.subplot(grid[-1, 0])
        ax2.matshow(    vls.reshape(1, len(vls)),
                        cmap=colors,
                        vmin=0,
                        vmax=vmax,
                        interpolation='bilinear',
                        )
        ax2.set_aspect(ax2.get_data_ratio() ** -1 * .5 / grid_wdth)
        ax2.tick_params(axis='x', bottom=True, top=False, labelbottom=True, labeltop=False, labelsize=label_size, direction='out')
        ax2.tick_params(axis='y', left=False, right=False, labelright=False, labelleft=False)
        ax2.set_xticks(tcks_loc)
        ax2.set_xticklabels(tcks)

        fig.get_axes()[0].set_title('Contacts freuqency')
        grid.tight_layout(fig)
        matplotlib.pyplot.savefig(fname + '.' + fmt, format=fmt)
        matplotlib.pyplot.close()

    def save_histo(self, fname, all_inds_stc2=True, fmt='svg'):
        """
        Saves histogram of contact counts for each atom.

        Arguments:
        fname -- str; name of file to be created.
        titles -- dict int: str; keys are indexes of histos, values are title to be set.
        all_inds_stc2 -- bool; True by default
        """
        max_bars = 15

        trg_vls_all = self.cmtx.sum(axis=1) / float(self.n)
        pep_vls = self.cmtx.sum(axis=0) / float(self.n)

        inds1 = numpy.nonzero(trg_vls_all)[0]
        if all_inds_stc2:
            inds2 = range(self.cmtx.shape[1])
        else:
            inds2 = numpy.nonzero(pep_vls)[0]

        trg_vls = list(numpy.array(trg_vls_all)[inds1,])
        trg_lbls = list(numpy.array(self.s1)[inds1,])

        pep_lbls = [self.s2[i] for i in inds2]

        max_y = max([.05] + trg_vls_all)

        chunks = _chunk_lst(trg_vls, max_bars, extend_last=0.)
        grid = matplotlib.pyplot.GridSpec(2 + len(chunks), 1)
        size = (10, 3 * len(chunks))
        fig = matplotlib.pyplot.figure(figsize=size)

        peptH = mk_histo(matplotlib.pyplot.subplot(grid[0, 0]), pep_vls, pep_lbls, ylim=(0, max([.05] + pep_vls)))[0]
        sbplts = [matplotlib.pyplot.subplot(grid[i, 0]) for i in range(1, len(chunks) + 1)]
        targBH = mk_histo(sbplts, chunks, _chunk_lst(trg_lbls, max_bars, extend_last=''), ylim=(0, max_y))
        targAH = mk_histo(matplotlib.pyplot.subplot(grid[-1, 0]), trg_vls_all, self.s1, ylim=(0, max_y))[0]

        peptH.set_title('Histogram of peptide contacts')
        targBH[0].set_title('Histogram of receptor contacts - detailed analysis')
        targAH.set_title('Histogram of receptor contacts - summary analysis')

        for sfig in targBH + [peptH, targAH]:
            sfig.set_ylabel('Contact frequency')
            sfig.set_xlabel('Residue id')

        xloc = targAH.get_xticks()
        if len(xloc) > max_bars:
            inds = numpy.linspace(0, len(xloc) - 1, max_bars).astype(int)
            targAH.set_xticks(xloc[inds,])
            targAH.set_xticklabels(numpy.array(self.s1)[inds,])

        grid.tight_layout(fig)
        matplotlib.pyplot.savefig(fname + '.' + fmt, format=fmt)
        matplotlib.pyplot.close()

    def save_txt(self, stream):
        """Saves contact list in CSV format.

        Argument:
        stream -- file-like object; stream to which text will be passed.
        """
        inds1, inds2 = numpy.nonzero(self.cmtx)
        stream.write("# n=%i\n" % self.n)
        for m1, m2, (c1, c2) in zip([self.s1[i] for i in inds1], [self.s2[i] for i in inds2], zip(inds1, inds2)):
           stream.write("%s\t%s\t%.3f\n" % (m1, m2, self.cmtx[c1, c2]))

    def save_all(self, fname, norm_n=False, break_long_x=50):
        """Creates txt and png of given name."""
        with open(fname + '.txt', 'w') as f:
            self.save_txt(f)
        self.save_fig(fname, norm_n=norm_n, break_long_x=break_long_x)

    def __add__(self, other):
        """Addition of cmaps sums their matrices.

        Raises ValueError for cmaps of different particles.
        """
        if self.s1 != other.s1 or self.s2 != other.s2:
            raise ValueError("Cannot sum different particles' contact maps.")
        return ContactMap(self.cmtx + other.cmtx, self.s1, self.s2, self.n + other.n)
