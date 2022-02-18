# Time-scrambling method

The time-scrambling algorithm is a widely used method for estimating a weighted reference isotropic exposure level with the
which compares the distribution of observed events.

In order to produce a data map and a reference map with equatorial coordinates binned into a fine (HEALPix) grid. 

The data map contains the directions of events in equatorial coordinates (![\alpha,\delta](https://latex.codecogs.com/svg.latex?&space;\alpha,\delta)))
calculated from their respetive local coordinates (![\theta, \phi](https://latex.codecogs.com/svg.latex?&space;\theta,\phi)) and arrival time t. The distribution of events in this
The map is still highly anisotropic as it preserves the uneven exposure of the detector to different parts of the sky.

The anisotropy in the data map has two sources: one due to detector exposure effects as mentioned above, and a second from
astrophysical origin. These two sources of anisotropy can be required stable in their corresponding frames: the sidereal frame
of astrophysical source, and the detector (or local) frame of anisotropy by exposure. The time-scrambling algorithm
untangles these two sources, keeping the local coordinates of each event fixed while mixing the astrophysical coordinates.
To do this, each event in the data set is assigned a random time taken from the time distribution of all events in the data set.
data in order to track deficiencies in data acquisition.
For each real event in the data, 20 ``false'' events are generated in order to reduce statistical fluctuations. 
The false events are used to build the reference map. This process effectively shuffles the right ascension ![\alpha](https://latex.codecogs.com/svg.latex?&space;\alpha) of the event
keeping local coordinates and detector rate
no change from actual data. This procedure destroys the autocorrelation
of events in sidereal coordinates that may have existed in the data and is an estimate of what the distribution of events would have
in the case that there was no sidereal anisotropy.

The amount by which each event is shuffled in right ascension depends on the length of the time window 
![\Delta t](https://latex.codecogs.com/svg.latex?&space;\Delta&space;t) from which event times are sampled. The length of the time window can be several minutes to a maximum of 24 hours. The length of the window indicates the time scale over which the detector acceptance is estimated to be stable in local coordinates and the events are not mixed over the times at which the detector changes direction.
setting. If the environmental conditions change very slowly and the detector operation is quite stable,
a time window of maximum 24 hours can be used.

The estimation of the reference level determined by the time randomization method has some disadvantages. At the South Pole, the declination of
an event ![\delta](https://latex.codecogs.com/svg.latex?&space;\delta) has degeneracy with the zenith angle ![\theta](https://latex.codecogs.com/svg.latex?&space;\theta) in local coordinates, regardless of the arrival time of the event. Since the time of events is
randomized, these remain in the same declination band, which reduces the sensitivity of the anisotropic search for structures oriented mostly along declination parallels.
For example, a dipole anisotropy aligned with the Earth's axis of rotation would produce a distribution of events that would be
indistinguishable from a zenith angle effect caused by detector using time randomization method and therefore not detected.
