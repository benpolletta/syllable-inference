classdef Sequence

    properties

        Words cell
        Sylbs cell
        SylbList cell
        Phones cell
        PhoneList cell
        PhoneDurations double
        SylbDurations double
        WordDurations double
        SpeechRate double

    end

    properties (Access = private)

        WordBoundaries cell
        SylbBoundaries cell
        PhonesPerWord cell
        PhonesPerSylb cell

    end

    methods

        function obj = Sequence(Words, PhoneDurations)

            if nargin < 2, PhoneDurations = []; end

            obj.Words = Words;
            obj.Sylbs = cellfun(@(x) split(x, '*'), obj.Words, 'UniformOutput', false);
            obj.SylbList = cat(1, obj.Sylbs{:});
            obj.Phones = cellfun(@(x) cellfun(@(y) split(y, '/'), x, 'UniformOutput', false), obj.Sylbs, 'UniformOutput', false);
            SylbPhoneList = cellfun(@(x) cat(1, x{:}), obj.Phones, 'unif', 0);
            obj.PhoneList = cat(1, SylbPhoneList{:});

            obj.PhonesPerSylb = cellfun(@(x) cellfun(@length, x), obj.Phones, 'UniformOutput', false);
            obj.PhonesPerWord = cellfun(@sum, obj.PhonesPerSylb, 'UniformOutput', false);

            obj.SylbBoundaries = counts2BoundsCellRecursive(obj.PhonesPerSylb);
            SylbBoundariesMat = cell2mat(obj.SylbBoundaries);
            obj.SylbBoundaries = mat2cell(SylbBoundariesMat, 2, ones(1, length(SylbBoundariesMat)));
            obj.WordBoundaries = counts2BoundsCellRecursive(obj.PhonesPerWord);

            if isempty(PhoneDurations)
    
                [~, num_phones] = cumsumCellRecursive(obj.PhonesPerWord);

                obj.PhoneDurations = nan(1, num_phones);

            end

            % This code calculates durations of syllables as a cell of
            % vectors, w/ syllables in the same word grouped in the same
            % cell.
            % obj.SylbDurations = cellfun(@(x)arrayfun(@(i) sum(obj.PhoneDurations(x(1, i):x(2, i))), 1:size(x, 2)), obj.SylbBoundaries, 'UniformOutput', false);
            % This code calculates durations of syllables as a vector.
            obj.SylbDurations = cellfun(@(x) sum(obj.PhoneDurations(x(1):x(2))), obj.SylbBoundaries);
            % obj.Sylb_durations = cellfun(@cell2mat, obj.Sylb_durations, 'UniformOutput', 0);
            % obj.Sylb_durations = cell2mat(obj.Sylb_durations);
            obj.WordDurations = cellfun(@(x) sum(obj.PhoneDurations(x(1):x(2))), obj.WordBoundaries);
            % obj.WordDurations = cell2mat(obj.WordDurations);
            obj.SpeechRate = mean(obj.WordDurations);

        end

        function obj = setRate()


        end

%         function bounds = countsToBounds(counts)
% 
%             cumulative_sum = cumsum(counts);
% 
%             lower_bounds = [0 cumulative_sum(1:(end - 1))] + 1;
%                 
%             bounds = [lower_bounds; cumulative_sum];
% 
%         end

    end

    methods (Access = private)

    end


end