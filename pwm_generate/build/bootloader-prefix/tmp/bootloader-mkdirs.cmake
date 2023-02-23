# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "C:/Users/15148/esp/esp-idf/components/bootloader/subproject"
  "D:/Github/CapSpeaker/pwm_generate/build/bootloader"
  "D:/Github/CapSpeaker/pwm_generate/build/bootloader-prefix"
  "D:/Github/CapSpeaker/pwm_generate/build/bootloader-prefix/tmp"
  "D:/Github/CapSpeaker/pwm_generate/build/bootloader-prefix/src/bootloader-stamp"
  "D:/Github/CapSpeaker/pwm_generate/build/bootloader-prefix/src"
  "D:/Github/CapSpeaker/pwm_generate/build/bootloader-prefix/src/bootloader-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "D:/Github/CapSpeaker/pwm_generate/build/bootloader-prefix/src/bootloader-stamp/${subDir}")
endforeach()
