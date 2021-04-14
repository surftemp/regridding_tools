# -*- coding: utf-8 -*-

#    sst-services
#    Copyright (C) 2020  National Centre for Earth Observation (NCEO)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

from .extractor import Extractor
from datetime import datetime, timedelta

class ExtractTimeSeriesTests(unittest.TestCase):

    def test_time_daily(self):
        expected1 = [
            (datetime(1999, 12, 30, 12, 0), datetime(1999, 12, 30, 12, 0), datetime(1999, 12, 30, 12, 0)),
            (datetime(1999, 12, 31, 12, 0), datetime(1999, 12, 31, 12, 0), datetime(1999, 12, 31, 12, 0))]

        expected2 = [
            (datetime(2000, 1, 1, 12, 0), datetime(2000, 1, 1, 12, 0), datetime(2000, 1, 1, 12, 0)),
            (datetime(2000, 1, 2, 12, 0), datetime(2000, 1, 2, 12, 0), datetime(2000, 1, 2, 12, 0))]

        periods1 = Extractor.createTimePeriods("daily",datetime(1999,12,30,12,0,0),datetime(1999,12,31,12,0,0))
        self.assertEqual(expected1,periods1)

        periods2 = Extractor.createTimePeriods("1", datetime(2000,1,1,12,0,0),datetime(2000,1,2,12,0,0))
        self.assertEqual(expected2, periods2)

    def test_time_pentad(self):
        expected = [
            (datetime(2003, 1, 1, 12, 0), datetime(2003, 1, 3, 12, 0), datetime(2003, 1, 5, 12, 0))]

        periods = Extractor.createTimePeriods("5-day", datetime(2003, 1, 1, 12, 0, 0),
                                                        datetime(2003, 1, 5, 12, 0, 0))
        self.assertEqual(expected, periods)

        # TODO add more tests

    def test_time_dekad(self):
        expected = [
            (datetime(2003, 1, 1, 12, 0), datetime(2003, 1, 5, 12, 0), datetime(2003, 1, 10, 12, 0))]

        periods = Extractor.createTimePeriods("10-day", datetime(2003, 1, 1, 12, 0, 0),
                                                        datetime(2003, 1, 10, 12, 0, 0))
        self.assertEqual(expected, periods)

        # TODO add more tests

    def test_time_monthly(self):
        expected = [
            (datetime(2003, 1, 1, 12, 0), datetime(2003, 1, 15, 12, 0), datetime(2003, 1, 31, 12, 0))]

        periods = Extractor.createTimePeriods("monthly", datetime(2003, 1, 1, 12, 0, 0),
                                                        datetime(2003, 1, 31, 12, 0, 0))
        self.assertEqual(expected, periods)

        # TODO add tests for february

    def test_time_ndaily_7(self):
        expected = [
            (datetime(2003, 1, 1, 12, 0), datetime(2003, 1, 4, 12, 0), datetime(2003, 1, 7, 12, 0)),
            (datetime(2003, 1, 8, 12, 0), datetime(2003, 1, 11, 12, 0), datetime(2003, 1, 14, 12, 0))
        ]

        periods = Extractor.createTimePeriods("7", datetime(2003, 1, 1, 12, 0, 0),
                                                        datetime(2003, 1, 14, 12, 0, 0))
        self.assertEqual(expected, periods)

    def test_time_ndaily_1(self):
        expected = [
            (datetime(2003, 1, 1, 12, 0)+timedelta(days=d),
             datetime(2003, 1, 1, 12, 0)+timedelta(days=d),
             datetime(2003, 1, 1, 12, 0)+timedelta(days=d)) for d in range(0,365)
        ]

        periods = Extractor.createTimePeriods("1", datetime(2003, 1, 1, 12, 0, 0),
                                                        datetime(2003, 12, 31, 12, 0, 0))
        self.assertEqual(expected, periods)

    def test_time_ndaily_2(self):

        periods1999 = Extractor.createTimePeriods("2", datetime(1999, 1, 1, 12, 0, 0),
                                                        datetime(1999, 12, 31, 12, 0, 0))
        expected_last1999 = (datetime(1999, 12, 31, 12, 0), datetime(1999, 12, 31, 12, 0), datetime(1999, 12, 31, 12, 0))

        self.assertEqual(periods1999[-1],expected_last1999)
        self.assertEqual(183,len(periods1999))

        periods2000 = Extractor.createTimePeriods("2", datetime(2000, 1, 1, 12, 0, 0),
                                                            datetime(2000, 12, 31, 12, 0, 0))
        expected_last2000 = (datetime(2000, 12, 30, 12, 0), datetime(2000, 12, 31, 0, 0),
                             datetime(2000, 12, 31, 12, 0))

        self.assertEqual(periods2000[-1], expected_last2000)
        self.assertEqual(183, len(periods2000))

    def test_time_ndaily_3(self):

        periods1999 = Extractor.createTimePeriods("3", datetime(1999, 1, 1, 12, 0, 0),
                                                        datetime(1999, 12, 31, 12, 0, 0))
        expected_last1999 = (datetime(1999, 12, 30, 12, 0), datetime(1999, 12, 31, 0, 0), datetime(1999, 12, 31, 12, 0))

        self.assertEqual(periods1999[-1],expected_last1999)
        self.assertEqual(122,len(periods1999))

        periods2000 = Extractor.createTimePeriods("3", datetime(2000, 1, 1, 12, 0, 0),
                                                            datetime(2000, 12, 31, 12, 0, 0))
        expected_last2000 = (datetime(2000, 12, 29, 12, 0), datetime(2000, 12, 30, 12, 0),
                             datetime(2000, 12, 31, 12, 0))

        self.assertEqual(periods2000[-1], expected_last2000)
        self.assertEqual(122, len(periods2000))

    def test_time_ndaily_4(self):

        periods1999 = Extractor.createTimePeriods("4", datetime(1999, 1, 1, 12, 0, 0),
                                                        datetime(1999, 12, 31, 12, 0, 0))
        expected_last1999 = (datetime(1999, 12, 27, 12, 0), datetime(1999, 12, 29, 12, 0), datetime(1999, 12, 31, 12, 0))

        self.assertEqual(periods1999[-1],expected_last1999)
        self.assertEqual(91,len(periods1999))

        periods2000 = Extractor.createTimePeriods("4", datetime(2000, 1, 1, 12, 0, 0),
                                                            datetime(2000, 12, 31, 12, 0, 0))
        expected_last2000 = (datetime(2000, 12, 30, 12, 0), datetime(2000, 12, 31, 0, 0),
                             datetime(2000, 12, 31, 12, 0))

        self.assertEqual(periods2000[-1], expected_last2000)
        self.assertEqual(92, len(periods2000))

    def test_time_ndaily_5(self):

        periods1999 = Extractor.createTimePeriods("5", datetime(1999, 1, 1, 12, 0, 0),
                                                        datetime(1999, 12, 31, 12, 0, 0))
        expected_last1999 = (datetime(1999, 12, 27, 12, 0), datetime(1999, 12, 29, 12, 0), datetime(1999, 12, 31, 12, 0))

        self.assertEqual(periods1999[-1],expected_last1999)
        self.assertEqual(73,len(periods1999))

        periods2000 = Extractor.createTimePeriods("5", datetime(2000, 1, 1, 12, 0, 0),
                                                            datetime(2000, 12, 31, 12, 0, 0))
        expected_last2000 = (datetime(2000, 12, 26, 12, 0), datetime(2000, 12, 29, 0, 0),
                             datetime(2000, 12, 31, 12, 0))

        self.assertEqual(periods2000[-1], expected_last2000)
        self.assertEqual(73, len(periods2000))

    def test_time_ndaily_6(self):

        periods1999 = Extractor.createTimePeriods("6", datetime(1999, 1, 1, 12, 0, 0),
                                                        datetime(1999, 12, 31, 12, 0, 0))
        expected_last1999 = (datetime(1999, 12, 27, 12, 0), datetime(1999, 12, 29, 12, 0), datetime(1999, 12, 31, 12, 0))

        self.assertEqual(periods1999[-1],expected_last1999)
        self.assertEqual(61,len(periods1999))

        periods2000 = Extractor.createTimePeriods("6", datetime(2000, 1, 1, 12, 0, 0),
                                                            datetime(2000, 12, 31, 12, 0, 0))
        expected_last2000 = (datetime(2000, 12, 26, 12, 0), datetime(2000, 12, 29, 0, 0),
                             datetime(2000, 12, 31, 12, 0))

        self.assertEqual(periods2000[-1], expected_last2000)
        self.assertEqual(61, len(periods2000))

    def test_time_ndaily_7(self):

        periods1999 = Extractor.createTimePeriods("7", datetime(1999, 1, 1, 12, 0, 0),
                                                        datetime(1999, 12, 31, 12, 0, 0))
        expected_last1999 = (datetime(1999, 12, 24, 12, 0), datetime(1999, 12, 28, 0, 0), datetime(1999, 12, 31, 12, 0))

        self.assertEqual(periods1999[-1],expected_last1999)
        self.assertEqual(52,len(periods1999))

        periods2000 = Extractor.createTimePeriods("7", datetime(2000, 1, 1, 12, 0, 0),
                                                            datetime(2000, 12, 31, 12, 0, 0))
        expected_last2000 = (datetime(2000, 12, 23, 12, 0), datetime(2000, 12, 27, 12, 0),
                             datetime(2000, 12, 31, 12, 0))

        self.assertEqual(periods2000[-1], expected_last2000)
        self.assertEqual(52, len(periods2000))

if __name__ == '__main__':
    unittest.main()