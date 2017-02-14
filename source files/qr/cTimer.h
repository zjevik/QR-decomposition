/**
 * \class arg::cTimer
 * \brief An utility timer.
 *
 * Can measure time range in seconds and millis.
 *
 * \author Pavel Kromer, (c) 2007 - 2011
 *
 * \version 2.0
 *
 * \b History:
 *		- initial version, 2007, pkromer
 * 		- doxygen comments, 26-07-2007, pkromer
 * 		- 2011-01, 	pkromer,	for ms precision used sys/time.h. The use of new methods Cpu* highly suggested.
 * 		- 2011-02,	pkromer,	old methods removed, made portable
 *
 */
#ifndef __CTIMER__
#define __CTIMER__

#ifdef __GNUC__
	#include <sys/time.h>
#endif

#include <ctime>

namespace arg
{
	class cTimer
	{
		#ifndef __GNUC__
			typedef struct timeval
			{
				unsigned long tv_sec;
				unsigned long tv_usec;
			};

			void gettimeofday(timeval * my_time, void * ignore) const
			{
				ignore = ignore;			// to get rid of unreferenced formal parameters
				(*my_time).tv_sec = (*my_time).tv_usec = 0;
				(*my_time).tv_sec = (unsigned long) time(NULL);
			}
		#endif

		protected:
			timeval m_TickerStart;
			unsigned long m_TickerIntervalMS;
			timeval m_TStart, m_TStop;

		public:
			inline void CpuStart(void); ///< Start interval on host.
			inline cTimer & CpuStop(void); ///< Stop interval on host.
			inline double CpuMillis(void); ///< Computes time interval on host.
			inline double CpuSeconds(void); ///< Computes time interval on host.

			inline void StartTicking(const unsigned long sec);
			inline bool IsTicking(void) const;
	};

	inline double cTimer::CpuMillis(void)
	{
		return ((m_TStop.tv_sec * 1e6 + m_TStop.tv_usec) - (m_TStart.tv_sec * 1e6 + m_TStart.tv_usec)) / 1e3;
	}

	inline double cTimer::CpuSeconds(void)
	{
		return CpuMillis() / 1e3;
	}

	inline void cTimer::CpuStart(void)
	{
		gettimeofday(&m_TStart, NULL);
	}

	inline cTimer& cTimer::CpuStop(void)
	{
		gettimeofday(&m_TStop, NULL);
		return *this;
	}

	inline void cTimer::StartTicking(const unsigned long sec)
	{
		gettimeofday(&m_TickerStart, NULL);
		m_TickerIntervalMS = sec * 1000;
	}

	inline bool cTimer::IsTicking(void) const
	{
		timeval now;
		gettimeofday(&now, NULL);

		double millis = ((now.tv_sec * 1e6 + now.tv_usec) - (m_TickerStart.tv_sec * 1e6 + m_TickerStart.tv_usec)) / 1e3;

		return (m_TickerIntervalMS > 0) && (millis < m_TickerIntervalMS);
	}
}

#endif
