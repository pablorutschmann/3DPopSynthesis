import Synthesis.post as post

TEST = post.population('SynthesisRuns/initsynth', 9)

TEST.final_mass_distribution(thresholds=[0,6.0e-11])
TEST.final_mass_distribution(cumulative=True, thresholds=[0,1.0e-11,9.0e-11])
TEST.final_mass_distribution(cumulative=False, density=True, thresholds=[1.0e-11,9.0e-11])
TEST.final_mass_distribution(cumulative=True, density=True, thresholds=[1.0e-11,9.0e-11])

TEST.final_radius_distribution(thresholds=[0,1.0e-11])
TEST.final_radius_distribution(cumulative=True, thresholds=[0,1.0e-11])
TEST.final_radius_distribution(cumulative=False, density=True, thresholds=[1.0e-11,9.0e-11])
TEST.final_radius_distribution(cumulative=True, density=True, thresholds=[1.0e-11,9.0e-11])

TEST.final_wmf_radial_distribution(thresholds=[0,1.0e-11])
TEST.final_wmf_radial_distribution(cumulative=True, thresholds=[0,1.0e-11])
TEST.final_wmf_radial_distribution(cumulative=False, density=True, thresholds=[1.0e-11,9.0e-11])
TEST.final_wmf_radial_distribution(cumulative=True, density=True, thresholds=[1.0e-11,9.0e-11])

