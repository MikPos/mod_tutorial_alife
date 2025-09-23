import collections
from typing import Dict, List, Tuple, Any, Optional
from dataclasses import dataclass, field
from mod import causality
import re
import time as t
import bisect
import copy
import numpy as np
from scipy import stats
from scipy.spatial.distance import cosine
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from sklearn.preprocessing import normalize
import pickle
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict


@dataclass
class ReactionStatistics:
    """
    Statistics for a single reaction across multiple simulation runs.
    
    Attributes:
        mean: Mean number of occurrences across simulations
        std: Standard deviation of occurrences
        min: Minimum number of occurrences
        max: Maximum number of occurrences
        median: Median number of occurrences
        raw_counts: List of raw counts from each simulation
        n_simulations: Number of simulations analyzed
    """
    mean: float
    std: float
    min: int
    max: int
    median: float
    raw_counts: List[int] = field(default_factory=list)
    n_simulations: int = 0


@dataclass
class TimePointStatistics:
    """
    Statistical measures for a species at a specific time point.
    
    Attributes:
        mean: Mean value across simulations
        std: Standard deviation
        sem: Standard error of the mean
        min: Minimum value
        max: Maximum value
        median: Median value
        q25: 25th percentile
        q75: 75th percentile
        ci_lower: Lower bound of confidence interval
        ci_upper: Upper bound of confidence interval
        n_samples: Number of samples
        raw_values: List of raw values from each simulation
    """
    mean: float
    std: float
    sem: float
    min: float
    max: float
    median: float
    q25: float
    q75: float
    ci_lower: float
    ci_upper: float
    n_samples: int
    raw_values: List[float] = field(default_factory=list)


@dataclass
class SpeciesTrajectoryStatistics:
    """
    Statistical summary for a species trajectory across multiple simulations.
    
    Attributes:
        times: List of time points
        species: Species name
        mean: List of mean values at each time point
        std: List of standard deviations at each time point
        sem: List of standard errors at each time point
        min: List of minimum values at each time point
        max: List of maximum values at each time point
        median: List of median values at each time point
        q25: List of 25th percentiles at each time point
        q75: List of 75th percentiles at each time point
        ci_lower: List of lower confidence interval bounds
        ci_upper: List of upper confidence interval bounds
        n_simulations: Number of simulations
        confidence_level: Confidence level used for intervals
    """
    times: List[float]
    species: str
    mean: List[float]
    std: List[float]
    sem: List[float]
    min: List[float]
    max: List[float]
    median: List[float]
    q25: List[float]
    q75: List[float]
    ci_lower: List[float]
    ci_upper: List[float]
    n_simulations: int
    confidence_level: float


@dataclass
class SummaryStatistics:
    """
    Summary statistics for a numerical field across multiple simulations.
    
    Attributes:
        mean: Mean value across simulations
        std: Standard deviation
        min: Minimum value
        max: Maximum value
        median: Median value
        raw_values: List of raw values from each simulation
    """
    mean: float
    std: float
    min: float
    max: float
    median: float
    raw_values: List[float] = field(default_factory=list)


@dataclass
class ReactionInfo:
    """
    Information about a reaction event.
    
    Attributes:
        sources: List of source molecule SMILES
        targets: List of target molecule SMILES
        rules: List of reaction rule names
    """
    sources: List[str]
    targets: List[str]
    rules: List[str]


@dataclass
class StateInfo:
    """
    Information about a simulation state.
    
    Attributes:
        time: Current simulation time
        iteration: Current iteration number
        time_increment: Time increment since last state
        event_type: Type of event that triggered the state change
        marking: Current marking as dictionary
        reaction: Optional reaction information if this was a reaction event
    """
    time: float
    iteration: int
    time_increment: float
    event_type: Optional[str]
    marking: Dict[str, int]
    reaction: Optional[ReactionInfo] = None


class SimulationCache:
    """
    A class to store and manage simulation results with convenience methods.
    
    This class provides structured storage for simulation data and methods
    to query and analyze the stored information.
    """
    
    def __init__(self):
        """Initialize the simulation cache with empty data structures."""
        self.times = []
        self.recomputes = []
        self.recomputes_avoided = []
        self.deadlocks = []
        self.states = []
        self.smiles = {}  # species name -> smiles
        self.parsed_markings = []  # Parsed version of markings as dictionaries
        self.reaction_types = collections.defaultdict(list)
        self.reaction_times = collections.defaultdict(list)
        self.parameters = {}
    
    def clear(self):
        """
        Clear all cached simulation data.
        
        Should be called before starting a new simulation if reusing the same cache.
        """
        return Stoch.__init__()
	
    def copy(self):
        """
        Create a deep copy of the simulation cache.
        
        Returns:
            SimulationCache: A deep copy of this cache instance
        """
        return copy.deepcopy(self)
    
    def get_time_stepped_markings(self, interval: float, start_time: Optional[float] = None, 
                                  end_time: Optional[float] = None) -> List[Tuple[float, Dict[str, int]]]:
        """
        Get markings at regular time intervals.
        
        For each time step, returns the last marking that occurred before or at that time.
        
        Args:
            interval: Time interval between steps
            start_time: Starting time (defaults to first time in cache)
            end_time: Ending time (defaults to last time in cache)
            
        Returns:
            List of (time, marking) tuples at regular intervals
            
        Examples:
            >>> cache = SimulationCache()
            >>> cache.times = [0.0, 0.5, 1.2, 2.1, 3.0]
            >>> cache.parsed_markings = [{'A': 1}, {'A': 2}, {'A': 1}, {'A': 3}, {'A': 0}]
            >>> stepped = cache.get_time_stepped_markings(1.0)
            >>> len(stepped)
            4
            >>> stepped[0]
            (0.0, {'A': 1})
            >>> stepped[1]
            (1.0, {'A': 1})
        """
        if not self.times or not self.parsed_markings:
            return []
        
        if start_time is None:
            start_time = self.times[0]
        if end_time is None:
            end_time = self.times[-1]
        
        result = []
        current_time = start_time
        
        while current_time <= end_time:
            # Find the index of the last time <= current_time
            idx = bisect.bisect_right(self.times, current_time) - 1
            
            if idx >= 0:
                marking = self.parsed_markings[idx]
                result.append((current_time, marking.copy()))
            
            current_time += interval
        
        return result
    
    def get_species_trajectory(self, species: str, interval: Optional[float] = None,
                               start_time: Optional[float] = None, 
                               end_time: Optional[float] = None) -> List[Tuple[float, int]]:
        """
        Get the trajectory of a specific species over time.
        
        Args:
            species: Species name (SMILES string or original name)
            interval: If provided, return values at regular intervals; otherwise return all recorded values
            start_time: Starting time (defaults to first time in cache)
            end_time: Ending time (defaults to last time in cache)
            
        Returns:
            List of (time, count) tuples for the specified species
            
        Examples:
            >>> cache = SimulationCache()
            >>> cache.times = [0.0, 1.0, 2.0]
            >>> cache.parsed_markings = [{'A': 5}, {'A': 3}, {'A': 1}]
            >>> traj = cache.get_species_trajectory('A')
            >>> traj
            [(0.0, 5), (1.0, 3), (2.0, 1)]
        """
        if interval is not None:
            stepped_markings = self.get_time_stepped_markings(interval, start_time, end_time)
            return [(time, marking.get(species, 0)) for time, marking in stepped_markings]
        else:
            if start_time is None:
                start_time = self.times[0] if self.times else 0.0
            if end_time is None:
                end_time = self.times[-1] if self.times else 0.0
            
            result = []
            for i, time in enumerate(self.times):
                if start_time <= time <= end_time:
                    count = self.parsed_markings[i].get(species, 0)
                    result.append((time, count))
            
            return result
    
    def get_all_species(self) -> set:
        """
        Get all species that appear in the simulation.
        
        Returns:
            Set of all species names (SMILES strings or original names)
            
        Examples:
            >>> cache = SimulationCache()
            >>> cache.parsed_markings = [{'A': 1, 'B': 2}, {'A': 0, 'C': 1}]
            >>> sorted(cache.get_all_species())
            ['A', 'B', 'C']
        """
        all_species = set()
        for marking in self.parsed_markings:
            all_species.update(marking.keys())
        return all_species
    
    def get_reaction_counts(self) -> Dict[str, int]:
        """
        Get the total number of times each reaction occurred.
        
        Returns:
            Dictionary mapping reaction names to occurrence counts
            
        Examples:
            >>> cache = SimulationCache()
            >>> cache.reaction_types = {'rxn1': [(0.1, 1), (0.2, 2)], 'rxn2': [(0.15, 1)]}
            >>> counts = cache.get_reaction_counts()
            >>> counts['rxn1']
            2
            >>> counts['rxn2']
            1
        """
        return {reaction: len(occurrences) for reaction, occurrences in self.reaction_types.items()}
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the simulation results.
        
        Returns:
            Dictionary containing simulation summary statistics
        """
        summary = {
            'total_iterations': len(self.times),
            'simulation_time': self.times[-1] - self.times[0] if len(self.times) > 1 else 0.0,
            'total_recomputes': len(self.recomputes),
            'total_recomputes_avoided': len(self.recomputes_avoided),
            'total_deadlocks': len(self.deadlocks),
            'unique_species': len(self.get_all_species()),
            'unique_reactions': len(self.reaction_types),
            'total_reactions': sum(len(occurrences) for occurrences in self.reaction_types.values())
        }
        
        if self.parameters:
            summary['parameters'] = self.parameters.copy()
        
        return summary

    def get_times(self):
        """Get the list of simulation times."""
        return self.times
	
    def get_recomputes(self):
        """Get the list of recompute iterations."""
        return self.recomputes
	
    def get_recomputes_avoided(self):
        """Get the list of avoided recompute iterations."""
        return self.recomputes_avoided
	
    def get_deadlocks(self):
        """Get the list of deadlock occurrences."""
        return self.deadlocks
	
    def get_states(self):
        """Get the list of simulation states."""
        return self.states
	
    def get_parsed_markings(self):
        """Get the list of parsed markings."""
        return self.parsed_markings
	
    def get_reaction_types(self):
        """Get the dictionary of reaction types and their occurrences."""
        return self.reaction_types


class SimulationStatistics:
    """
    A class to compute statistics across multiple simulation cache instances.
    
    This class provides methods to analyze multiple simulation runs and compute
    statistics like mean, standard deviation, confidence intervals, etc.
    """
    
    def __init__(self, caches: List[SimulationCache]):
        """
        Initialize with a list of simulation cache instances.
        
        Args:
            caches: List of SimulationCache instances to analyze
        """
        self.caches = caches
        if not caches:
            raise ValueError("At least one cache instance is required")
    
    def get_time_stepped_statistics(self, interval: float, start_time: Optional[float] = None, 
                                   end_time: Optional[float] = None, 
                                   confidence_level: float = 0.95) -> Dict[str, Any]:
        """
        Compute statistics for species concentrations at regular time intervals across all simulations.
        
        Args:
            interval: Time interval between steps
            start_time: Starting time (defaults to earliest start across all caches)
            end_time: Ending time (defaults to latest end across all caches)
            confidence_level: Confidence level for confidence intervals (default 0.95)
            
        Returns:
            Dictionary containing time points and statistics for each species
            
        Examples:
            >>> cache1 = SimulationCache()
            >>> cache1.times = [0.0, 1.0, 2.0]
            >>> cache1.parsed_markings = [{'A': 10}, {'A': 8}, {'A': 6}]
            >>> cache2 = SimulationCache()
            >>> cache2.times = [0.0, 1.0, 2.0]
            >>> cache2.parsed_markings = [{'A': 12}, {'A': 10}, {'A': 8}]
            >>> stats = SimulationStatistics([cache1, cache2])
            >>> result = stats.get_time_stepped_statistics(1.0)
            >>> 'times' in result
            True
            >>> 'species_stats' in result
            True
        """
        # Determine time bounds across all caches
        if start_time is None:
            start_time = min(cache.times[0] for cache in self.caches if cache.times)
        if end_time is None:
            end_time = max(cache.times[-1] for cache in self.caches if cache.times)
        
        # Get all species across all caches
        all_species = set()
        for cache in self.caches:
            all_species.update(cache.get_all_species())
        
        # Get time-stepped markings for each cache
        all_stepped_markings = []
        for cache in self.caches:
            stepped = cache.get_time_stepped_markings(interval, start_time, end_time)
            all_stepped_markings.append(stepped)
        
        if not all_stepped_markings or not all_stepped_markings[0]:
            return {'times': [], 'species_stats': {}}
        
        # Extract time points (should be the same for all caches)
        time_points = [t for t, _ in all_stepped_markings[0]]
        
        # Organize data by species and time point
        species_stats = {}
        for species in all_species:
            species_data = []
            
            # For each time point, collect values across all simulations
            for t_idx in range(len(time_points)):
                values = []
                for cache_markings in all_stepped_markings:
                    if t_idx < len(cache_markings):
                        _, marking = cache_markings[t_idx]
                        values.append(marking.get(species, 0))
                    else:
                        values.append(0)  # Default to 0 if time point not available
                species_data.append(values)
            
            # Compute statistics for each time point
            time_stats = []
            for values in species_data:
                if values:
                    values_array = np.array(values)
                    mean_val = np.mean(values_array)
                    std_val = np.std(values_array, ddof=1) if len(values_array) > 1 else 0.0
                    
                    # Compute confidence interval
                    if len(values_array) > 1:
                        sem = stats.sem(values_array)  # Standard error of the mean
                        h = sem * stats.t.ppf((1 + confidence_level) / 2., len(values_array) - 1)
                        ci_lower = mean_val - h
                        ci_upper = mean_val + h
                    else:
                        ci_lower = ci_upper = mean_val
                    
                    time_stats.append(TimePointStatistics(
                        mean=mean_val,
                        std=std_val,
                        sem=stats.sem(values_array) if len(values_array) > 1 else 0.0,
                        min=np.min(values_array),
                        max=np.max(values_array),
                        median=np.median(values_array),
                        q25=np.percentile(values_array, 25),
                        q75=np.percentile(values_array, 75),
                        ci_lower=ci_lower,
                        ci_upper=ci_upper,
                        n_samples=len(values_array),
                        raw_values=values
                    ))
                else:
                    # No data available
                    time_stats.append(TimePointStatistics(
                        mean=0.0, std=0.0, sem=0.0,
                        min=0.0, max=0.0, median=0.0,
                        q25=0.0, q75=0.0,
                        ci_lower=0.0, ci_upper=0.0,
                        n_samples=0, raw_values=[]
                    ))
            
            species_stats[species] = time_stats
        
        return {
            'times': time_points,
            'species_stats': species_stats,
            'n_simulations': len(self.caches),
            'confidence_level': confidence_level,
            'interval': interval
        }
    
    def get_species_trajectory_statistics(self, species: str, interval: float,
                                        start_time: Optional[float] = None,
                                        end_time: Optional[float] = None,
                                        confidence_level: float = 0.95) -> SpeciesTrajectoryStatistics:
        """
        Get statistical summary for a specific species trajectory across all simulations.
        
        Args:
            species: Species name to analyze
            interval: Time interval between steps
            start_time: Starting time
            end_time: Ending time
            confidence_level: Confidence level for confidence intervals
            
        Returns:
            SpeciesTrajectoryStatistics object with times and statistical measures for the species
        """
        stats_data = self.get_time_stepped_statistics(interval, start_time, end_time, confidence_level)
        
        if species not in stats_data['species_stats']:
            return SpeciesTrajectoryStatistics(
                times=stats_data['times'],
                species=species,
                mean=[],
                std=[],
                sem=[],
                min=[],
                max=[],
                median=[],
                q25=[],
                q75=[],
                ci_lower=[],
                ci_upper=[],
                n_simulations=stats_data['n_simulations'],
                confidence_level=confidence_level
            )
        
        species_data = stats_data['species_stats'][species]
        
        return SpeciesTrajectoryStatistics(
            times=stats_data['times'],
            species=species,
            mean=[point.mean for point in species_data],
            std=[point.std for point in species_data],
            sem=[point.sem for point in species_data],
            min=[point.min for point in species_data],
            max=[point.max for point in species_data],
            median=[point.median for point in species_data],
            q25=[point.q25 for point in species_data],
            q75=[point.q75 for point in species_data],
            ci_lower=[point.ci_lower for point in species_data],
            ci_upper=[point.ci_upper for point in species_data],
            n_simulations=stats_data['n_simulations'],
            confidence_level=confidence_level
        )
    
    def get_reaction_statistics(self) -> Dict[str, ReactionStatistics]:
        """
        Get statistics for reaction counts across all simulations.
        
        Returns:
            Dictionary mapping reaction names to ReactionStatistics objects
        """
        # Get all unique reactions
        all_reactions = set()
        for cache in self.caches:
            all_reactions.update(cache.reaction_types.keys())
        
        reaction_stats = {}
        for reaction in all_reactions:
            counts = []
            for cache in self.caches:
                count = len(cache.reaction_types.get(reaction, []))
                counts.append(count)
            
            if counts:
                counts_array = np.array(counts)
                reaction_stats[reaction] = ReactionStatistics(
                    mean=np.mean(counts_array),
                    std=np.std(counts_array, ddof=1) if len(counts_array) > 1 else 0.0,
                    min=int(np.min(counts_array)),
                    max=int(np.max(counts_array)),
                    median=np.median(counts_array),
                    raw_counts=counts,
                    n_simulations=len(counts)
                )
        
        return reaction_stats
    
    def get_summary_statistics(self) -> Dict[str, Any]:
        """
        Get summary statistics across all simulation runs.
        
        Returns:
            Dictionary with comprehensive summary statistics. The 'summary_stats' key contains
            SummaryStatistics objects for numerical fields, while other keys contain metadata.
        """
        summaries = [cache.get_summary() for cache in self.caches]
        
        # Aggregate numerical values
        numerical_fields = ['total_iterations', 'simulation_time', 'total_recomputes', 
                          'total_recomputes_avoided', 'total_deadlocks', 'unique_species',
                          'unique_reactions', 'total_reactions']
        
        aggregated = {}
        for field in numerical_fields:
            values = [s[field] for s in summaries if field in s]
            if values:
                values_array = np.array(values)
                aggregated[field] = SummaryStatistics(
                    mean=np.mean(values_array),
                    std=np.std(values_array, ddof=1) if len(values_array) > 1 else 0.0,
                    min=np.min(values_array),
                    max=np.max(values_array),
                    median=np.median(values_array),
                    raw_values=values
                )
        
        return {
            'n_simulations': len(self.caches),
            'summary_stats': aggregated,
            'individual_summaries': summaries
        }
    
    def plot_concentration_profiles(self, 
                                   plot_species: Optional[List[str]] = None,
                                   output_path: Optional[str] = None,
                                   interval: float = 1.0,
                                   start_time: Optional[float] = None,
                                   end_time: Optional[float] = None,
                                   uncertainty_type: str = 'std',
                                   confidence_level: float = 0.95,
                                   zoom_factor: float = 0.2,
                                   padding_factor: float = 1.05,
                                   figsize: Tuple[float, float] = (14, 8)) -> None:
        """
        Plot concentration profiles of selected species with uncertainty bands.
        
        Creates a dual-panel plot showing concentration profiles with a main plot
        (zoomed to show detail) and an unzoomed view. Supports both standard deviation
        and confidence intervals for uncertainty visualization.
        
        Args:
            plot_species: List of species names to plot. If None, all species will be plotted.
            output_path: Path to save the plot (should include file extension). If None, the plot will be saved to the current working directory.
            interval: Time interval between data points (default 1.0)
            start_time: Starting time for analysis (defaults to earliest across all caches)
            end_time: Ending time for analysis (defaults to latest across all caches)
            uncertainty_type: Type of uncertainty to plot ('std' for standard deviation, 
                            'ci' for confidence intervals)
            confidence_level: Confidence level for confidence intervals (default 0.95)
            zoom_factor: Factor to zoom the main plot to visualize species with lower concentrations. (default 0.2 = 20% of max)
            padding_factor: Factor to pad the y-axis, this is used to ensure that all the data is visible. (default 1.05)
            figsize: Figure size as (width, height) tuple (default (14, 8))
            
        Raises:
            ValueError: If uncertainty_type is not 'std' or 'ci', or if no species to plot
            FileNotFoundError: If output directory doesn't exist
            
        Examples:
            >>> # Create simulation statistics from caches
            >>> stats = SimulationStatistics([cache1, cache2, cache3])
            >>> 
            >>> # Plot with standard deviation uncertainty
            >>> stats.plot_concentration_profiles(
            ...     plot_species=['Glycolaldehyde', 'Glyoxylate'],
            ...     output_path='concentration_profiles.png',
            ...     uncertainty_type='std'
            ... )
            >>> 
            >>> # Plot with confidence intervals
            >>> stats.plot_concentration_profiles(
            ...     plot_species=['Glycolaldehyde', 'Glyoxylate'],
            ...     output_path='concentration_profiles_ci.png',
            ...     uncertainty_type='ci',
            ...     confidence_level=0.99
            ... )
        """

        if uncertainty_type not in ['std', 'ci']:
            raise ValueError("uncertainty_type must be 'std' or 'ci'")
        
        # Ensure output directory exists
        if output_path:
            output_dir = os.path.dirname(output_path)
        else:
            output_dir = os.getcwd()
            output_path = os.path.join(output_dir, 'concentration_profiles.png')
        if output_dir and not os.path.exists(output_dir):
            raise FileNotFoundError(f"Output directory does not exist: {output_dir}")
        
        # Get time-stepped statistics
        stats_data = self.get_time_stepped_statistics(
            interval=interval,
            start_time=start_time,
            end_time=end_time,
            confidence_level=confidence_level
        )
        
        times = np.array(stats_data['times'])
        species_stats = stats_data['species_stats']
        
        # Filter to only include requested species that exist in the data
        if plot_species is None:
            plot_species = list(species_stats.keys())
        available_species = [s for s in plot_species if s in species_stats]
        if not available_species:
            raise ValueError(f"None of the requested species {plot_species} found in simulation data")
        
        if len(available_species) != len(plot_species):
            missing = set(plot_species) - set(available_species)
            print(f"Warning: Species not found in data: {missing}")
        
        # Create figure with specified size
        fig = plt.figure(figsize=figsize)
        
        # Create main plot (left side, larger)
        ax1 = plt.axes([0.1, 0.1, 0.6, 0.8])  # [left, bottom, width, height]
        
        # Create unzoomed plot (right side, smaller)
        ax2 = plt.axes([0.75, 0.1, 0.25, 0.48])  # [left, bottom, width, height]
        
        # Create colormap with distinct colors
        colors = plt.cm.tab10(np.linspace(0, 1, len(available_species)))
        
        # Calculate global maximum for y-axis scaling (including uncertainty bounds)
        global_max = 0
        for species in available_species:
            for i in range(len(species_stats[species])):
                point = species_stats[species][i]
                if uncertainty_type == 'std':
                    # For standard deviation: mean + std
                    max_value = point.mean + point.std
                else:  # confidence intervals
                    # For confidence intervals: use the upper bound
                    max_value = point.ci_upper
                global_max = max(global_max, max_value)
        max_time = times[-1] if len(times) > 0 else 1.0
        
        # Get SMILES dictionary from first cache
        smiles_dict = {}
        if self.caches:
            smiles_dict = self.caches[0].smiles
        
        # Plot each species
        for i, species in enumerate(available_species):
            # Extract data for this species
            means = np.array([point.mean for point in species_stats[species]])
            
            # Calculate uncertainty based on type
            if uncertainty_type == 'std':
                uncertainties = np.array([point.std for point in species_stats[species]])
            else:  # confidence intervals
                uncertainties_lower = np.array([point.ci_lower for point in species_stats[species]])
                uncertainties_upper = np.array([point.ci_upper for point in species_stats[species]])
            
            # Create label with SMILES if available
            if smiles_dict and species in smiles_dict and smiles_dict[species]:
                label = f"{smiles_dict[species]}"
            else:
                label = species
            
            # Plot on main axis (zoomed)
            ax1.plot(times, means, label=label, color=colors[i], linewidth=1.5)
            if uncertainty_type == 'std':
                ax1.fill_between(times, means - uncertainties, means + uncertainties, 
                               color=colors[i], alpha=0.2)
            else:
                ax1.fill_between(times, uncertainties_lower, uncertainties_upper, 
                               color=colors[i], alpha=0.2)
            
            # Plot on unzoomed axis
            ax2.plot(times, means, color=colors[i], linewidth=1.5)
            if uncertainty_type == 'std':
                ax2.fill_between(times, means - uncertainties, means + uncertainties, 
                               color=colors[i], alpha=0.2)
            else:
                ax2.fill_between(times, uncertainties_lower, uncertainties_upper, 
                               color=colors[i], alpha=0.2)
        
        # Configure main plot (zoomed)
        uncertainty_label = 'Standard Deviation' if uncertainty_type == 'std' else f'{confidence_level*100:.0f}% Confidence Interval'
        ax1.set_title(f'Molecule Count Profiles ({zoom_factor*100:.0f}% of max)\nUncertainty: {uncertainty_label}', 
                     fontsize=15)
        ax1.set_xlabel('Simulation Time', fontsize=15)
        ax1.set_ylabel('Average Molecule Count', fontsize=15)
        ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, global_max * zoom_factor * padding_factor)
        ax1.set_xlim(0, max_time)
        
        # Configure unzoomed plot
        ax2.set_title('Full Range View', fontsize=15)
        ax2.set_xlabel('Simulation Time', fontsize=15)
        ax2.set_ylabel('Average Molecule Count', fontsize=15)
        ax2.set_ylim(0, global_max*padding_factor)
        ax2.set_xlim(0, max_time)
        ax2.grid(True, alpha=0.3)
        ax2.tick_params(labelsize=8)  # Make labels smaller for compact view
        
        # Save the plot
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Concentration profiles plot saved to: {output_path}")
        print(f"Plotted {len(available_species)} species with {uncertainty_label.lower()}")
    
    def calculate_species_similarity(self, 
                                   interval: float = 1.0,
                                   start_time: Optional[float] = None,
                                   end_time: Optional[float] = None) -> pd.DataFrame:
        """
        Calculate cosine similarities between species based on their mean concentration profiles.
        
        Args:
            interval: Time interval between data points (default 1.0)
            start_time: Starting time for analysis (defaults to earliest across all caches)
            end_time: Ending time for analysis (defaults to latest across all caches)
            
        Returns:
            pd.DataFrame: DataFrame with cosine similarity scores for each species pair
            
        Examples:
            >>> stats = SimulationStatistics([cache1, cache2, cache3])
            >>> similarity_matrix = stats.calculate_species_similarity()
            >>> print(similarity_matrix.shape)
            (n_species, n_species)
        """
        # Get time-stepped statistics
        stats_data = self.get_time_stepped_statistics(
            interval=interval,
            start_time=start_time,
            end_time=end_time
        )
        
        species_stats = stats_data['species_stats']
        species_list = list(species_stats.keys())
        
        if len(species_list) < 2:
            raise ValueError("Need at least 2 species to calculate similarities")
        
        # Filter out species with all-zero time series
        non_zero_species = []
        for species in species_list:
            means = np.array([point.mean for point in species_stats[species]])
            if np.any(means > 0):  # At least one non-zero value
                non_zero_species.append(species)
        
        if len(non_zero_species) < 2:
            raise ValueError("Need at least 2 species with non-zero concentrations to calculate similarities")
        
        if len(non_zero_species) < len(species_list):
            zero_species = set(species_list) - set(non_zero_species)
            print(f"Warning: Excluding {len(zero_species)} species with all-zero concentrations: {list(zero_species)}")
        
        # Create matrix of mean concentration profiles for non-zero species only
        matrix = np.vstack([np.array([point.mean for point in species_stats[species]]) 
                           for species in non_zero_species]).T
        
        # Normalize the matrix so each species vector has unit norm
        matrix_norm = normalize(matrix, axis=0)
        
        # Calculate cosine similarity matrix
        n_species = len(non_zero_species)
        sim_matrix = np.zeros((n_species, n_species))
        
        for i in range(n_species):
            for j in range(n_species):
                # Cosine similarity = 1 - cosine distance
                sim_matrix[i, j] = 1 - cosine(matrix_norm[:, i], matrix_norm[:, j])
        
        # Create a DataFrame from the similarity matrix
        sim_df = pd.DataFrame(sim_matrix, index=non_zero_species, columns=non_zero_species)
        
        return sim_df
    
    def cluster_species_by_similarity(self, 
                                    similarity_matrix: Optional[pd.DataFrame] = None,
                                    output_path: Optional[str] = None,
                                    n_clusters: int = 3,
                                    interval: float = 1.0,
                                    start_time: Optional[float] = None,
                                    end_time: Optional[float] = None) -> Dict[str, int]:
        """
        Cluster species based on pairwise similarity using hierarchical clustering.
        
        Args:
            similarity_matrix: Pre-computed similarity matrix (optional)
            output_path: Path to save dendrogram plot (optional)
            n_clusters: Target number of clusters (default 3)
            interval: Time interval for similarity calculation (default 1.0)
            start_time: Starting time for analysis (defaults to earliest across all caches)
            end_time: Ending time for analysis (defaults to latest across all caches)
            
        Returns:
            Dict[str, int]: Mapping from species name to cluster ID (1-based)
            
        Examples:
            >>> stats = SimulationStatistics([cache1, cache2, cache3])
            >>> clusters = stats.cluster_species_by_similarity(n_clusters=4)
            >>> print(f"Species A is in cluster {clusters['A']}")
        """
        # Calculate similarity matrix if not provided
        if similarity_matrix is None:
            similarity_matrix = self.calculate_species_similarity(
                interval=interval,
                start_time=start_time,
                end_time=end_time
            )
        
        # Convert similarities to distances
        distances = 1 - similarity_matrix.values
        
        # Convert to condensed distance matrix for linkage
        condensed_distances = squareform(distances, checks=False)
        
        # Perform hierarchical clustering
        linkage_matrix = linkage(condensed_distances, method='ward')
        
        # Extract cluster assignments
        cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
        
        # Create species to cluster mapping
        species_list = list(similarity_matrix.index)
        clusters = {species: int(label) for species, label in zip(species_list, cluster_labels)}
        
        # Plot dendrogram if output path provided
        if output_path:
            plt.figure(figsize=(12, 8))
            dendrogram(linkage_matrix, labels=species_list, leaf_rotation=90)
            plt.title(f'Species Clustering Dendrogram (n_clusters={n_clusters})')
            plt.xlabel('Species')
            plt.ylabel('Distance')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Dendrogram saved to: {output_path}")
        
        return clusters
    
    def plot_species_time_heatmap(self, 
                                 clusters: Optional[Dict[str, int]] = None,
                                 output_path: Optional[str] = None,
                                 top_n: int = 50,
                                 interval: float = 1.0,
                                 start_time: Optional[float] = None,
                                 end_time: Optional[float] = None,
                                 figsize: Tuple[float, float] = (14, 10)) -> None:
        """
        Create a heatmap of normalized species concentrations over time, with species
        ordered by clustering results or abundance.
        
        Args:
            clusters: Dictionary mapping species names to cluster IDs (optional)
            output_path: Path to save the heatmap (optional)
            top_n: Number of most abundant species to include (default 50)
            interval: Time interval between data points (default 1.0)
            start_time: Starting time for analysis (defaults to earliest across all caches)
            end_time: Ending time for analysis (defaults to latest across all caches)
            figsize: Figure size as (width, height) tuple (default (14, 10))
            
        Examples:
            >>> stats = SimulationStatistics([cache1, cache2, cache3])
            >>> clusters = stats.cluster_species_by_similarity(n_clusters=4)
            >>> stats.plot_species_time_heatmap(clusters=clusters, top_n=20)
        """
        # Get time-stepped statistics
        stats_data = self.get_time_stepped_statistics(
            interval=interval,
            start_time=start_time,
            end_time=end_time
        )
        
        times = np.array(stats_data['times'])
        species_stats = stats_data['species_stats']
        
        if not species_stats:
            print("No species data available for heatmap.")
            return
        
        # Get SMILES dictionary from first cache
        smiles_dict = {}
        if self.caches:
            smiles_dict = self.caches[0].smiles
        
        # Calculate maximum mean abundance for each species
        max_mean_abundance = {}
        for species, stats_list in species_stats.items():
            max_mean_abundance[species] = max([point.mean for point in stats_list])
        
        # Filter out species with zero abundance and sort by maximum mean abundance
        non_zero_species = [s for s in species_stats.keys() if max_mean_abundance[s] > 0]
        top_species = sorted(non_zero_species, 
                           key=lambda s: max_mean_abundance[s], reverse=True)[:top_n]
        
        if len(non_zero_species) < len(species_stats):
            zero_species = set(species_stats.keys()) - set(non_zero_species)
            print(f"Note: Excluding {len(zero_species)} species with zero concentrations from heatmap")
        
        print(f"Creating heatmap with top {len(top_species)} most abundant species")
        
        # Create heatmap data matrix
        heatmap_data = []
        species_labels = []
        
        # Check if we have any non-zero data
        has_data = False
        
        for species in top_species:
            means = np.array([point.mean for point in species_stats[species]])
            max_val = np.max(means)
            
            if max_val > 0:
                has_data = True
                # Normalize by maximum value for this species
                normalized_means = means / max_val
                heatmap_data.append(normalized_means)
                
                # Create label with SMILES if available
                if smiles_dict and species in smiles_dict and smiles_dict[species]:
                    label = f"{species}\n{smiles_dict[species]}"
                else:
                    label = species
                species_labels.append(label)
        
        if not has_data:
            print("No non-zero data available for heatmap.")
            return
        
        # Convert to numpy array
        heatmap_array = np.array(heatmap_data)
        
        # Create the heatmap
        plt.figure(figsize=figsize)
        
        # Create heatmap
        im = plt.imshow(heatmap_array, aspect='auto', cmap='viridis', interpolation='nearest')
        
        # Set labels and title
        plt.xlabel('Time Points', fontsize=12)
        plt.ylabel('Species', fontsize=12)
        plt.title(f'Normalized Species Concentration Heatmap (Top {len(top_species)} Species)', fontsize=14)
        
        # Set y-axis labels
        plt.yticks(range(len(species_labels)), species_labels, fontsize=8)
        
        # Set x-axis labels (show every nth time point to avoid crowding)
        n_time_points = len(times)
        if n_time_points > 20:
            step = max(1, n_time_points // 10)
            x_ticks = range(0, n_time_points, step)
            x_labels = [f'{times[i]:.1f}' for i in x_ticks]
        else:
            x_ticks = range(n_time_points)
            x_labels = [f'{times[i]:.1f}' for i in x_ticks]
        
        plt.xticks(x_ticks, x_labels, fontsize=10)
        
        # Add colorbar
        cbar = plt.colorbar(im)
        cbar.set_label('Normalized Concentration', fontsize=12)
        
        # Add cluster information if available
        if clusters:
            # Group species by cluster for better organization
            cluster_groups = defaultdict(list)
            for i, species in enumerate(top_species):
                if species in clusters:
                    cluster_id = clusters[species]
                    cluster_groups[cluster_id].append(i)
            
            # # Add cluster separators
            # y_positions = []
            # for cluster_id in sorted(cluster_groups.keys()):
            #     positions = cluster_groups[cluster_id]
            #     if positions:
            #         y_positions.append((min(positions), max(positions)))
            
            # # Draw cluster separators
            # for i, (start, end) in enumerate(y_positions[:-1]):
            #     plt.axhline(y=end + 0.5, color='red', linestyle='--', alpha=0.7, linewidth=1)
        
        plt.tight_layout()
        
        # Save plot if output path provided
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Species time heatmap saved to: {output_path}")
        else:
            plt.show()
        
        plt.close()
    
    def statistical_analysis(self, 
                              output_dir: str,
                              plot_species: Optional[List[str]] = None,
                              n_clusters: int = 3,
                              top_n_heatmap: int = 50,
                              interval: float = 1.0,
                              start_time: Optional[float] = None,
                              end_time: Optional[float] = None,
                              uncertainty_type: str = 'std',
                              confidence_level: float = 0.95,
                              zoom_factor: float = 0.2,
                              padding_factor: float = 1.05) -> Dict[str, Any]:
        """
        Perform a statistical analysis including similarity calculation, clustering, and plotting.
        
        This method combines multiple analysis functions to provide a complete overview
        of the simulation results, including:
        - Species concentration profiles with uncertainty
        - Species similarity matrix and clustering
        - Species time heatmap
        - Summary statistics
        
        Args:
            output_dir: Directory to save all plots and results
            plot_species: List of species to plot in concentration profiles (optional)
            n_clusters: Number of clusters for species clustering (default 3)
            top_n_heatmap: Number of top species for heatmap (default 50)
            interval: Time interval between data points (default 1.0)
            start_time: Starting time for analysis (defaults to earliest across all caches)
            end_time: Ending time for analysis (defaults to latest across all caches)
            uncertainty_type: Type of uncertainty ('std' or 'ci')
            confidence_level: Confidence level for confidence intervals (default 0.95)
            zoom_factor: Zoom factor for concentration profiles (default 0.2)
            padding_factor: Padding factor for y-axis (default 1.05)
            
        Returns:
            Dict containing analysis results including similarity matrix, clusters, and summary
            
        Examples:
            >>> stats = SimulationStatistics([cache1, cache2, cache3])
            >>> results = stats.comprehensive_analysis(
            ...     output_dir='analysis_results',
            ...     plot_species=['Glycolaldehyde', 'Glyoxylate'],
            ...     n_clusters=4
            ... )
            >>> print(f"Found {len(results['clusters'])} clusters")
        """
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"Starting comprehensive analysis...")
        print(f"Output directory: {output_dir}")
        print(f"Number of simulations: {len(self.caches)}")
        
        results = {}
        
        try:
            # 1. Calculate species similarity matrix
            print("\n1. Calculating species similarity matrix...")
            similarity_matrix = self.calculate_species_similarity(
                interval=interval,
                start_time=start_time,
                end_time=end_time
            )
            results['similarity_matrix'] = similarity_matrix
            print(f"   Computed similarity matrix for {len(similarity_matrix)} species")
            # 2. Cluster species by similarity
            print("\n2. Clustering species by similarity...")
            clusters = self.cluster_species_by_similarity(
                similarity_matrix=similarity_matrix,
                output_path=os.path.join(output_dir, 'species_clusters_dendrogram.png'),
                n_clusters=n_clusters
            )
            results['clusters'] = clusters
            print(f"   Clustered {len(clusters)} species into {n_clusters} clusters")
            
            # 3. Plot species time heatmap
            print("\n3. Creating species time heatmap...")
            self.plot_species_time_heatmap(
                clusters=clusters,
                output_path=os.path.join(output_dir, 'species_time_heatmap.png'),
                top_n=top_n_heatmap,
                interval=interval,
                start_time=start_time,
                end_time=end_time
            )
            
            # 4. Plot concentration profiles
            print("\n4. Plotting concentration profiles...")
            self.plot_concentration_profiles(
                plot_species=plot_species,
                output_path=os.path.join(output_dir, 'concentration_profiles.png'),
                interval=interval,
                start_time=start_time,
                end_time=end_time,
                uncertainty_type=uncertainty_type,
                confidence_level=confidence_level,
                zoom_factor=zoom_factor,
                padding_factor=padding_factor
            )
            
            # 5. Get summary statistics
            print("\n5. Computing summary statistics...")
            summary_stats = self.get_summary_statistics()
            results['summary_stats'] = summary_stats
            print(f"   Computed summary statistics for {summary_stats['n_simulations']} simulations")
            
            # 6. Save results to files
            print("\n6. Saving results...")
            
            # Save similarity matrix
            similarity_matrix.to_csv(os.path.join(output_dir, 'species_similarity_matrix.csv'))
            
            # Save clusters with SMILES
            clusters_df = pd.DataFrame(list(clusters.items()), columns=['species', 'cluster'])
            if self.caches and self.caches[0].smiles:
                clusters_df['smiles'] = clusters_df['species'].map(
                    lambda s: self.caches[0].smiles.get(s, "")
                )
            clusters_df = clusters_df.sort_values(['cluster', 'species'])
            clusters_df.to_csv(os.path.join(output_dir, 'species_clusters.csv'), index=False)
            
            # Save summary statistics
            summary_data = {
                'n_simulations': summary_stats['n_simulations'],
                'summary_stats': {k: {
                    'mean': v.mean,
                    'std': v.std,
                    'min': v.min,
                    'max': v.max,
                    'median': v.median
                } for k, v in summary_stats['summary_stats'].items()}
            }
            
            print(f"  Results saved to: {output_dir}")
            print(f"  Files created:")
            print(f"    - species_similarity_matrix.csv")
            print(f"    - species_clusters.csv")
            print(f"    - species_clusters_dendrogram.png")
            print(f"    - species_time_heatmap.png")
            print(f"    - concentration_profiles.png")
            print(f"    - summary_statistics.json")
            
        except Exception as e:
            print(f"\n Error during statistical analysis: {e}")
            import traceback
            traceback.print_exc()
            raise
        
        return results




def clear_cache():
    """
    Clear the global simulation cache.
    
    Should be called before starting a new simulation if reusing the same process.
    This function is provided for backward compatibility.
    
    Examples:
        >>> from callbacks import clear_cache, simulation_cache
        >>> simulation_cache.times.append(1.0)
        >>> len(simulation_cache.times)
        1
        >>> clear_cache()
        >>> len(simulation_cache.times)
        0
    """
    global simulation_cache
    simulation_cache.clear()


def parse_marking(marking) -> Dict[str, int]:
    """
    Parse the string representation of a Marking object into a dictionary.
    Uses direct string slicing for robustness.
    
    Args:
        marking: Marking object or its string representation
        
    Returns:
        Dict[str, int]: Dictionary mapping molecule names to counts
        
    Examples:
        >>> parse_marking("Marking{Glycolaldehyde: 9, Glyoxylate: 99, Mol-3: 1}")
        {'Glycolaldehyde': 9, 'Glyoxylate': 99, 'Mol-3': 1}
    """
    # Convert marking to string if it's not already
    marking_str = str(marking)
    
    # Check if the string starts with "Marking{"
    if not marking_str.startswith("Marking{"):
        return {}
    
    # Extract the content between the braces
    content = marking_str[8:-1]  # Remove "Marking{" and "}"
    
    # If the content is empty, return an empty dictionary
    if not content:
        return {}
    
    # Split the content by commas and process each key-value pair
    result = {}
    pairs = []
    
    # Handle commas inside molecule names with curly braces or other special characters
    # Start parsing character by character
    current_pair = ""
    brace_level = 0
    for char in content:
        if char == '{':
            brace_level += 1
            current_pair += char
        elif char == '}':
            brace_level -= 1
            current_pair += char
        elif char == ',' and brace_level == 0:
            # Only split on commas that are not inside braces
            pairs.append(current_pair.strip())
            current_pair = ""
        else:
            current_pair += char
    
    # Add the last pair if it exists
    if current_pair:
        pairs.append(current_pair.strip())
    
    # Process each key-value pair
    for pair in pairs:
        # Find the position of the colon
        colon_pos = pair.find(':')
        if colon_pos != -1:
            key = pair[:colon_pos].strip()
            value_str = pair[colon_pos+1:].strip()
            try:
                value = int(value_str)
                result[key] = value
            except ValueError:
                # Skip if the value is not an integer
                pass
    
    return result

# Global cache instance for backward compatibility
simulation_cache = SimulationCache()

def setCallbacks(sim, params=None, verbose=False, cache=None):
    """
    Set up callbacks for the stochastic simulator.
    
    Args:
        sim: A stochastic simulator instance from mod.stochsim
        params: Optional dictionary of simulation parameters to store with the results
        verbose: If True, print detailed output during simulation
        cache: Optional SimulationCache instance to use (defaults to global cache)
        
    Returns:
        SimulationCache: The cache instance being used
    """
    global simulation_cache 
    if cache is None:
        cache = simulation_cache
    # Clear the cache
    simulation_cache = SimulationCache()
    
    # Store simulation parameters if provided
    if params:
        cache.parameters = params.copy()
    
    # Capture initial state at time 0
    initial_marking_str = str(sim._marking)
    initial_parsed_marking = parse_marking(initial_marking_str)
    
    # Store SMILES for initial molecules
    for marking in initial_parsed_marking:
        if marking not in cache.smiles:
            for vertex in sim._dg.vertices:
                if str(vertex.graph)[1:-1] == marking:
                    cache.smiles[marking] = vertex.graph.smiles
                    break
    
    # Create initial smiles marking
    initial_smiles_marking = {}
    for marking in initial_parsed_marking:
        if marking in cache.smiles:
            initial_smiles_marking[cache.smiles[marking]] = initial_parsed_marking[marking]
        else:
            initial_smiles_marking[marking] = initial_parsed_marking[marking]
    
    # Store initial state information
    initial_state_info = StateInfo(
        time=0.0,
        iteration=0,
        time_increment=0.0,
        event_type=None,
        marking=initial_smiles_marking
    )
    
    cache.states.append(initial_state_info)
    cache.parsed_markings.append(initial_smiles_marking)
    cache.times.append(0.0)
    
    if verbose:
        print("Initial marking at t=0:", initial_marking_str)
    
    def onRecompute(i):
        """
        Called when a recomputation is performed.
        
        Args:
            i: Iteration number
        """
        cache.recomputes.append(i)
        if verbose:
            print("Recompute:", i)
    
    def onRecomputeAvoided(i):
        """
        Called when a recomputation is avoided.
        
        Args:
            i: Iteration number
        """
        cache.recomputes_avoided.append(i)
        if verbose:
            print("RecomputeAvoided:", i)
    
    def onDeadlock(time, i, trace):
        """
        Called when a deadlock is detected.
        
        Args:
            time: Current simulation time
            i: Iteration number
            trace: Simulation trace
        """
        cache.deadlocks.append((time, i))
        if verbose:
            print("Deadlock: {}, t={}".format(i, time))
    
    def onNewState(time, i, trace, e, tInc):
        """
        Called when a new state is reached.
        
        Args:
            time: Current simulation time
            i: Iteration number
            trace: Simulation trace
            e: Event that triggered the state change
            tInc: Time increment since last state
            
        Returns:
            True to continue simulation, False to stop
        """
        # Cache state information
        state_info = StateInfo(
            time=time,
            iteration=i,
            time_increment=tInc,
            event_type=str(type(e).__name__) if e else None,
            marking={}  # Will be filled later
        )
        cache.states.append(state_info)

        # Store the marking and parse it
        marking_str = str(sim._marking)
        parsed_marking = parse_marking(marking_str)

        for marking in parsed_marking:
            if marking not in cache.smiles:
                for vertex in sim._dg.vertices:
                    if str(vertex.graph)[1:-1] == marking:
                        cache.smiles[marking] = vertex.graph.smiles
                        break

        smiles_marking = {}
        # Replace internal names with smiles
        for marking in parsed_marking:
            if marking in cache.smiles:
                smiles_marking[cache.smiles[marking]] = parsed_marking[marking]
            else:
                smiles_marking[marking] = parsed_marking[marking]

        cache.parsed_markings.append(smiles_marking)
        cache.times.append(time)
        
        # Add the parsed marking to the state info
        state_info.marking = smiles_marking
        
        if verbose:
            print("Marking:", marking_str)
        
        # If this is an edge action (reaction), cache reaction type information
        if e and type(e) is causality.EdgeAction:
            for rule in e.edge.rules:
                rule_name = str(rule)
                cache.reaction_types[rule_name].append((time, i))
                cache.reaction_times[rule_name].append(time)
                
            # Store information about the reaction
            state_info.reaction = ReactionInfo(
                sources=[s.graph.smiles for s in e.edge.sources],
                targets=[t.graph.smiles for t in e.edge.targets],
                rules=[str(r) for r in e.edge.rules]
            )
        
        if verbose:
            print("New state: {}, t={}, delta t={}".format(i, time, tInc))
            print("  e={}".format(e))
        
        return True
        
    # Assign callbacks to simulator
    sim.onRecompute = onRecompute
    sim.onRecomputeAvoided = onRecomputeAvoided
    sim.onDeadlock = onDeadlock
    sim.onNewState = onNewState
    
    return cache


def load_legacy_cache(filepath: str) -> SimulationCache:
    """
    Load a legacy pickle file (old dictionary format) and convert it to a SimulationCache instance.
    
    Args:
        filepath: Path to the legacy pickle file
        
    Returns:
        SimulationCache: A new cache instance with the loaded data
        
    Examples:
        >>> cache = load_legacy_cache('doe_cache/mu10.0_sig1.0_l06.0_n3.0_glycol200_glyox2000_runs10_cache.pkl')
        >>> len(cache.times) > 0
        True
        >>> isinstance(cache.parameters, dict)
        True
    """
    with open(filepath, 'rb') as f:
        legacy_data = pickle.load(f)
    
    # Create a new cache instance
    cache = SimulationCache()
    
    # Copy over the data, handling any missing fields gracefully
    cache.times = legacy_data.get('times', [])
    cache.recomputes = legacy_data.get('recomputes', [])
    cache.recomputes_avoided = legacy_data.get('recomputes_avoided', [])
    cache.deadlocks = legacy_data.get('deadlocks', [])
    cache.states = legacy_data.get('states', [])
    cache.smiles = legacy_data.get('smiles', {})
    cache.parsed_markings = legacy_data.get('parsed_markings', [])
    cache.parameters = legacy_data.get('parameters', {})
    
    # Handle reaction types and times (these might be defaultdicts in legacy)
    cache.reaction_types = collections.defaultdict(list)
    if 'reaction_types' in legacy_data:
        cache.reaction_types.update(legacy_data['reaction_types'])
    
    cache.reaction_times = collections.defaultdict(list)
    if 'reaction_times' in legacy_data:
        cache.reaction_times.update(legacy_data['reaction_times'])
    
    return cache


def save_cache(cache: SimulationCache, filepath: str):
    """
    Save a SimulationCache instance to a pickle file.
    
    Args:
        cache: SimulationCache instance to save
        filepath: Path where to save the cache
        
    Examples:
        >>> cache = SimulationCache()
        >>> cache.times = [0.0, 1.0, 2.0]
        >>> save_cache(cache, 'my_cache.pkl')
        >>> loaded_cache = load_cache('my_cache.pkl')
        >>> loaded_cache.times == cache.times
        True
    """
    # Ensure directory exists
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    
    with open(filepath, 'wb') as f:
        pickle.dump(cache, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_cache(filepath: str) -> SimulationCache:
    """
    Load a SimulationCache instance from a pickle file.
    
    Args:
        filepath: Path to the saved cache file
        
    Returns:
        SimulationCache: The loaded cache instance
        
    Examples:
        >>> cache = load_cache('my_cache.pkl')
        >>> isinstance(cache, SimulationCache)
        True
    """
    with open(filepath, 'rb') as f:
        return pickle.load(f)


def save_multiple_caches(caches: List[SimulationCache], filepath: str, 
                        metadata: Optional[Dict[str, Any]] = None):
    """
    Save multiple SimulationCache instances for statistical analysis.
    
    Args:
        caches: List of SimulationCache instances
        filepath: Path where to save the cache collection
        metadata: Optional metadata about the simulation runs
        
    Examples:
        >>> caches = [SimulationCache(), SimulationCache()]
        >>> save_multiple_caches(caches, 'multi_run_study.pkl', {'description': 'Parameter sweep'})
        >>> loaded_caches, meta = load_multiple_caches('multi_run_study.pkl')
        >>> len(loaded_caches) == 2
        True
    """
    # Ensure directory exists
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    
    data = {
        'caches': caches,
        'metadata': metadata or {},
        'n_runs': len(caches),
        'saved_at': t.time()
    }
    
    with open(filepath, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_multiple_caches(filepath: str) -> Tuple[List[SimulationCache], Dict[str, Any]]:
    """
    Load multiple SimulationCache instances saved together.
    
    Args:
        filepath: Path to the saved cache collection
        
    Returns:
        Tuple of (caches_list, metadata_dict)
        
    Examples:
        >>> caches, metadata = load_multiple_caches('multi_run_study.pkl')
        >>> isinstance(caches, list)
        True
        >>> isinstance(metadata, dict)
        True
    """
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
    
    # Handle both new format and simple list format
    if isinstance(data, dict) and 'caches' in data:
        return data['caches'], data.get('metadata', {})
    elif isinstance(data, list):
        # Simple list of caches
        return data, {}
    else:
        raise ValueError(f"Unexpected data format in file {filepath}")


def load_legacy_cache_directory(directory: str, pattern: str = "*.pkl") -> List[SimulationCache]:
    """
    Load all legacy cache files from a directory.
    
    Args:
        directory: Directory containing legacy pickle files
        pattern: File pattern to match (default: "*.pkl")
        
    Returns:
        List of SimulationCache instances loaded from the directory
        
    Examples:
        >>> caches = load_legacy_cache_directory('doe_cache')
        >>> len(caches) > 0
        True
        >>> all(isinstance(cache, SimulationCache) for cache in caches)
        True
    """
    directory_path = Path(directory)
    cache_files = list(directory_path.glob(pattern))
    
    caches = []
    for cache_file in cache_files:
        try:
            cache = load_legacy_cache(str(cache_file))
            caches.append(cache)
            print(f"Loaded legacy cache: {cache_file.name}")
        except Exception as e:
            print(f"Failed to load {cache_file.name}: {e}")
    
    return caches


def convert_legacy_directory(input_directory: str, output_directory: str, 
                           pattern: str = "*.pkl"):
    """
    Convert all legacy cache files in a directory to the new format.
    
    Args:
        input_directory: Directory containing legacy pickle files
        output_directory: Directory to save converted files
        pattern: File pattern to match (default: "*.pkl")
        
    Examples:
        >>> convert_legacy_directory('doe_cache', 'converted_cache')
        # Converts all .pkl files from doe_cache to new format in converted_cache
    """
    input_path = Path(input_directory)
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    
    cache_files = list(input_path.glob(pattern))
    
    for cache_file in cache_files:
        try:
            # Load legacy cache
            cache = load_legacy_cache(str(cache_file))
            
            # Save in new format
            output_file = output_path / cache_file.name
            save_cache(cache, str(output_file))
            
            print(f"Converted: {cache_file.name} -> {output_file}")
        except Exception as e:
            print(f"Failed to convert {cache_file.name}: {e}")


def create_statistics_from_legacy_directory(directory: str, pattern: str = "*.pkl") -> SimulationStatistics:
    """
    Create a SimulationStatistics instance from all legacy cache files in a directory.
    
    Args:
        directory: Directory containing legacy pickle files
        pattern: File pattern to match (default: "*.pkl")
        
    Returns:
        SimulationStatistics: Statistics instance ready for analysis
        
    Examples:
        >>> stats = create_statistics_from_legacy_directory('doe_cache')
        >>> isinstance(stats, SimulationStatistics)
        True
        >>> stats.caches  # List of loaded caches
    """
    caches = load_legacy_cache_directory(directory, pattern)
    if not caches:
        raise ValueError(f"No cache files found in directory {directory}")
    
    return SimulationStatistics(caches)
